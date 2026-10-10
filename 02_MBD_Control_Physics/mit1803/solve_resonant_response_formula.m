function [sol_analytic, sys_tf, sim_data] = solve_resonant_response_formula(m, b, k, F0, omega, x0, v0, t_span)
% SOLVE_RESONANT_RESPONSE_FORMULA Solves m*x'' + b*x' + k*x = F0*cos(omega*t).
% Automatically switches between ERF (non-resonant) and RRF (resonant response).
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, lsim, pole).
%
% % MBD Concept Mapping: Encapsulates operator p(D) into continuous LTI tf(num, den).
% %                      Extracts system poles and implements Resonant Response Formula (RRF)
% %                      to capture secular growth x_p = B*t*exp(rt)/p'(r) under pole excitation.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 15 Resonance Benchmark):
%   m      - Inertia / Mass coefficient (default: 1.0)
%   b      - Damping coefficient (default: 0.0)
%   k      - Spring constant (default: 4.0)
%   F0     - Forcing amplitude (default: 1.0)
%   omega  - Forcing circular frequency (default: 2.0 rad/s)
%   x0     - Initial position x(0) (default: 0.0)
%   v0     - Initial velocity x'(0) (default: 0.0)
%   t_span - [t_start, t_end] simulation interval (default: [0, 20.0])
%
% Outputs:
%   sol_analytic - Struct containing closed-form function handles and resonance status
%   sys_tf       - Transfer function object G(s) from Control System Toolbox
%   sim_data     - Struct containing lsim numerical verification results [t, y, u]

    arguments
        m (1,1) double {mustBePositive} = 1.0
        b (1,1) double {mustBeNonnegative} = 0.0
        k (1,1) double {mustBePositive} = 4.0
        F0 (1,1) double {mustBeNumeric} = 1.0
        omega (1,1) double {mustBePositive} = 2.0
        x0 (1,1) double {mustBeNumeric} = 0.0
        v0 (1,1) double {mustBeNumeric} = 0.0
        t_span (1,2) double {mustBeNumeric} = [0, 20.0]
    end

    %% 1. LTI System Construction via Control System Toolbox
    den = [m, b, k];
    sys_tf = tf(1, den);
    sys_poles = pole(sys_tf);

    %% 2. Resonance Check & Particular Solution (ERF vs. RRF)
    % Characteristic polynomial: p(s) = m*s^2 + b*s + k
    % p'(s) = 2*m*s + b
    s_eval = 1i * omega;
    p_val  = m * (s_eval^2) + b * s_eval + k;
    p_der  = 2 * m * s_eval + b;

    if abs(p_val) < 1e-9
        % Pure Resonance -> Apply Resonant Response Formula (RRF)
        sol_analytic.is_resonant = true;
        
        % z_p(t) = F0 * t * exp(i*omega*t) / p'(i*omega)
        % For m*x'' + k*x = F0*cos(omega*t) where b = 0, p'(i*omega) = 2*m*i*omega
        % z_p(t) = F0 * t * (cos(wt) + i*sin(wt)) / (2*m*i*omega) = (F0*t / 2*m*w) * (sin(wt) - i*cos(wt))
        % x_p(t) = Re(z_p) = (F0 / (2*m*omega)) * t * sin(omega*t)
        Z_scale = F0 / p_der; % Complex scalar for t*exp(1i*w*t)
        x_p_handle = @(t) real(Z_scale .* t .* exp(1i * omega * t));
        v_p_handle = @(t) real(Z_scale .* (1 + 1i * omega * t) .* exp(1i * omega * t));
        
        sol_analytic.secular_rate = abs(Z_scale);
        sol_analytic.envelope_handle = @(t) abs(Z_scale) .* t;
        
    else
        % Non-Resonant -> Apply Standard ERF
        sol_analytic.is_resonant = false;
        Z_ss = F0 / p_val;
        x_p_handle = @(t) real(Z_ss * exp(1i * omega * t));
        v_p_handle = @(t) real(1i * omega * Z_ss * exp(1i * omega * t));
        
        sol_analytic.secular_rate = 0.0;
        sol_analytic.envelope_handle = @(t) abs(Z_ss) * ones(size(t));
    end

    %% 3. Homogeneous Superposition to Match Initial Conditions
    x_p_0 = x_p_handle(0);
    v_p_0 = v_p_handle(0);
    x_h_0 = x0 - x_p_0;
    v_h_0 = v0 - v_p_0;

    disc = b^2 - 4 * m * k;
    if disc < -1e-12
        alpha = -b / (2 * m);
        w_d = sqrt(4 * m * k - b^2) / (2 * m);
        c1 = x_h_0;
        c2 = (v_h_0 - alpha * c1) / w_d;
        x_h_handle = @(t) exp(alpha * t) .* (c1 * cos(w_d * t) + c2 * sin(w_d * t));
    elseif abs(disc) <= 1e-12
        r = -b / (2 * m);
        c1 = x_h_0;
        c2 = v_h_0 - r * c1;
        x_h_handle = @(t) (c1 + c2 * t) .* exp(r * t);
    else
        r1 = (-b + sqrt(disc)) / (2 * m);
        r2 = (-b - sqrt(disc)) / (2 * m);
        c2 = (v_h_0 - r1 * x_h_0) / (r2 - r1);
        c1 = x_h_0 - c2;
        x_h_handle = @(t) c1 * exp(r1 * t) + c2 * exp(r2 * t);
    end

    x_total_handle = @(t) x_p_handle(t) + x_h_handle(t);

    sol_analytic.x_total = x_total_handle;
    sol_analytic.x_particular = x_p_handle;
    sol_analytic.x_homogeneous = x_h_handle;

    %% 4. Numerical Cross-Check via lsim (Control System Toolbox)
    t_eval = linspace(t_span(1), t_span(2), 1000);
    u_eval = F0 * cos(omega * t_eval);
    
    % Convert tf to ss to properly handle non-zero initial conditions in lsim
    sys_ss = ss(sys_tf);
    
    % State transformation: for transfer function 1/(m*s^2 + b*s + k), canonical state is [x; v]
    % Compute response using lsim
    if x0 == 0 && v0 == 0
        y_sim = lsim(sys_ss, u_eval, t_eval);
    else
        % Fallback for general initial state coordinates
        ode_eq = @(t, X) [X(2); (F0*cos(omega*t) - b*X(2) - k*X(1))/m];
        [~, X_eval] = ode45(ode_eq, t_eval, [x0; v0]);
        y_sim = X_eval(:, 1);
    end

    sim_data.t = t_eval;
    sim_data.y = y_sim;
    sim_data.u = u_eval;

    % Validation error assessment
    x_ana_eval = x_total_handle(t_eval);
    max_error = max(abs(x_ana_eval(:) - y_sim(:)));
    sol_analytic.verification_max_error = max_error;

    %% 5. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Resonant Response Solver: MIT 18.03 Lecture 15   \n');
    fprintf('====================================================\n');
    fprintf('Plant Model             : %.2f*x'''' + %.2f*x'' + %.2f*x = %.2f*cos(%.2f*t)\n', ...
            m, b, k, F0, omega);
    fprintf('Resonance Detected      : %d\n', sol_analytic.is_resonant);
    if sol_analytic.is_resonant
        fprintf('Solution Type           : Secular Growth via RRF: x_p = B*t*exp(rt)/p''(r)\n');
        fprintf('Secular Amplitude Rate  : %.4f per second\n', sol_analytic.secular_rate);
    else
        fprintf('Solution Type           : Bounded Periodic Steady-State via ERF\n');
    end
    fprintf('System Poles            : \n');
    disp(sys_poles);
    fprintf('Verification Max Error  : %e\n', max_error);
    fprintf('====================================================\n\n');

    %% 6. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 950, 500]);
        hold on; grid on; box on;

        plot(t_eval, x_total_handle(t_eval), 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Analytic Solution x(t)');
        plot(t_eval, y_sim, 'ro', 'MarkerIndices', 1:40:length(t_eval), ...
             'MarkerSize', 4, 'DisplayName', 'Control Toolbox lsim Check');

        if sol_analytic.is_resonant
            plot(t_eval, sol_analytic.envelope_handle(t_eval), 'r--', 'LineWidth', 1.4, ...
                 'DisplayName', 'Linear Envelope \pm (t / 4)');
            plot(t_eval, -sol_analytic.envelope_handle(t_eval), 'r--', 'LineWidth', 1.4, ...
                 'HandleVisibility', 'off');
        end

        xlabel('Time t (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('State Position x(t)', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('Forced Dynamic Response (Resonance = %d, \\omega = %.2f rad/s)', ...
              sol_analytic.is_resonant, omega), 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 9);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end