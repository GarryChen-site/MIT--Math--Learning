function [sol_analytic, sys_metrics] = solve_second_order_forced_response(m, b, k, F0, omega, theta_in, x0, v0, t_span)
% SOLVE_SECOND_ORDER_FORCED_RESPONSE Solves m*x'' + b*x' + k*x = F0*cos(omega*t + theta_in).
%
% % Required Toolboxes: 
% %   - Base MATLAB (Zero external dependency).
%
% % MBD Concept Mapping: Demonstrates frequency response operator G(jw) = 1/p(jw).
% %                      Separates steady-state forced tracking x_p(t) via ERF
% %                      from transient natural decay x_h(t) via Superposition II.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 13 Section [4] Benchmark):
%   m        - Mass / Inertia coefficient (default: 1.0)
%   b        - Damping coefficient (default: 1.0)
%   k        - Spring stiffness coefficient (default: 2.0)
%   F0       - Amplitude of sinusoidal forcing (default: 1.0)
%   omega    - Excitation angular frequency in rad/s (default: 1.0)
%   theta_in - Input phase angle in radians (default: 0.0)
%   x0       - Initial displacement state x(0) (default: 0.0)
%   v0       - Initial velocity state x'(0) (default: 0.0)
%   t_span   - [t_start, t_end] simulation interval (default: [0, 16.0])
%
% Outputs:
%   sol_analytic - Struct containing function handles:
%                  .x_total, .x_steady_state, .x_transient
%   sys_metrics  - Struct containing ERF evaluation, amplitude, phase lag, and error

    arguments
        m (1,1) double {mustBePositive} = 1.0
        b (1,1) double {mustBeNonnegative} = 1.0
        k (1,1) double {mustBePositive} = 2.0
        F0 (1,1) double {mustBeNumeric} = 1.0
        omega (1,1) double {mustBePositive} = 1.0
        theta_in (1,1) double {mustBeNumeric} = 0.0
        x0 (1,1) double {mustBeNumeric} = 0.0
        v0 (1,1) double {mustBeNumeric} = 0.0
        t_span (1,2) double {mustBeNumeric} = [0, 16.0]
    end

    %% 1. Particular Solution via Exponential Response Formula (ERF)
    % Characteristic polynomial: p(s) = m*s^2 + b*s + k
    % Evaluate at driving frequency s = 1i*omega:
    p_val = m * (1i * omega)^2 + b * (1i * omega) + k;
    
    if abs(p_val) < 1e-12
        error('Pure Undamped Resonance: p(j*omega) = 0. ERF denominator vanishes.');
    end
    
    % Input phasor in polar form: U_in = F0 * exp(1i * theta_in)
    U_in = F0 * exp(1i * theta_in);
    
    % Output complex phasor: Z_p = U_in / p(1i*omega)
    Z_p = U_in / p_val;
    
    A_ss = abs(Z_p);
    phi_lag = -angle(Z_p); % Phase delay relative to cos(omega*t)
    
    % Particular solution and its analytical derivative
    x_p_handle = @(t) real(Z_p * exp(1i * omega * t));
    v_p_handle = @(t) real(1i * omega * Z_p * exp(1i * omega * t));
    
    x_p_0 = x_p_handle(0);
    v_p_0 = v_p_handle(0);

    %% 2. Homogeneous Transient Solution (Superposition II)
    disc = b^2 - 4 * m * k;
    
    % Remaining initial states to be satisfied by the homogeneous solution
    x_h_0 = x0 - x_p_0;
    v_h_0 = v0 - v_p_0;
    
    if disc < -1e-12
        % Underdamped regime
        alpha = -b / (2 * m);
        w_d = sqrt(4 * m * k - b^2) / (2 * m);
        
        c1 = x_h_0;
        c2 = (v_h_0 - alpha * c1) / w_d;
        
        x_h_handle = @(t) exp(alpha * t) .* (c1 * cos(w_d * t) + c2 * sin(w_d * t));
        
    elseif abs(disc) <= 1e-12
        % Critically damped regime
        r = -b / (2 * m);
        c1 = x_h_0;
        c2 = v_h_0 - r * c1;
        
        x_h_handle = @(t) (c1 + c2 * t) .* exp(r * t);
        
    else
        % Overdamped regime
        r1 = (-b + sqrt(disc)) / (2 * m);
        r2 = (-b - sqrt(disc)) / (2 * m);
        
        c2 = (v_h_0 - r1 * x_h_0) / (r2 - r1);
        c1 = x_h_0 - c2;
        
        x_h_handle = @(t) c1 * exp(r1 * t) + c2 * exp(r2 * t);
    end

    % Total closed-form analytical response
    x_total_handle = @(t) x_p_handle(t) + x_h_handle(t);

    sol_analytic.x_total = x_total_handle;
    sol_analytic.x_steady_state = x_p_handle;
    sol_analytic.x_transient = x_h_handle;
    sol_analytic.Z_phasor = Z_p;

    %% 3. Metrics Packaging
    sys_metrics.w_n = sqrt(k / m);
    sys_metrics.zeta = b / (2 * sqrt(m * k));
    sys_metrics.p_eval = p_val;
    sys_metrics.steady_amplitude = A_ss;
    sys_metrics.steady_phase_lag_rad = phi_lag;
    sys_metrics.steady_phase_lag_deg = phi_lag * (180 / pi);
    sys_metrics.gain = A_ss / F0;

    %% 4. Numerical Cross-Check via ode45
    t_eval = linspace(t_span(1), t_span(2), 600);
    ode_fun = @(t, X) [X(2); (F0 * cos(omega * t + theta_in) - b * X(2) - k * X(1)) / m];
    
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [t_num, X_num] = ode45(ode_fun, t_eval, [x0; v0], opts);
    
    x_ana_eval = x_total_handle(t_num);
    max_error = max(abs(x_ana_eval - X_num(:, 1)));
    sys_metrics.verification_max_error = max_error;

    %% 5. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Second-Order Forced ODE Solver: MIT 18.03 Lec 13\n');
    fprintf('====================================================\n');
    fprintf('Model: %.2f*x'''' + %.2f*x'' + %.2f*x = %.2f*cos(%.2f*t + %.2f)\n', ...
            m, b, k, F0, omega, theta_in);
    fprintf('Natural Frequency (w_n) : %.4f rad/s\n', sys_metrics.w_n);
    fprintf('Damping Ratio (zeta)    : %.4f\n', sys_metrics.zeta);
    fprintf('ERF Characteristic p(jw): %+.4f %+.4fi\n', real(p_val), imag(p_val));
    fprintf('Steady-State Phasor Z_p : %+.4f %+.4fi\n', real(Z_p), imag(Z_p));
    fprintf('Steady Amplitude (A_ss) : %.4f\n', sys_metrics.steady_amplitude);
    fprintf('Steady Phase Lag (phi)  : %.4f rad (%.2f degrees)\n', ...
            sys_metrics.steady_phase_lag_rad, sys_metrics.steady_phase_lag_deg);
    fprintf('Verification Max Error  : %e\n', max_error);
    fprintf('====================================================\n\n');

    %% 6. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 950, 500]);
        hold on; grid on; box on;

        u_input = F0 * cos(omega * t_eval + theta_in);
        plot(t_eval, u_input, 'k:', 'LineWidth', 1.2, ...
             'DisplayName', sprintf('Forcing Input u(t) (Amp = %.1f)', F0));
        plot(t_eval, x_total_handle(t_eval), 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Total Response x(t)');
        plot(t_eval, x_p_handle(t_eval), 'r--', 'LineWidth', 1.6, ...
             'DisplayName', sprintf('Steady-State x_p(t) (Amp = %.3f, Lag = %.1f^\\circ)', A_ss, sys_metrics.steady_phase_lag_deg));
        plot(t_eval, x_h_handle(t_eval), 'g-.', 'LineWidth', 1.4, ...
             'DisplayName', 'Decaying Transient x_h(t)');
        plot(t_num(1:25:end), X_num(1:25:end, 1), 'mo', 'MarkerSize', 5, ...
             'DisplayName', 'ode45 Numerical Check');

        xlabel('Time t (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('State Position x(t)', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('Forced Dynamics: \\omega = %.2f rad/s, A_{ss} = %.3f, Lag = %.1f^\\circ', ...
              omega, A_ss, sys_metrics.steady_phase_lag_deg), 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end