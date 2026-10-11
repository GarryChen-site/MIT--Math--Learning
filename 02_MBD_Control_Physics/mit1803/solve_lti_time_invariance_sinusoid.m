function [sol_analytic, sys_tf, sim_data] = solve_lti_time_invariance_sinusoid(a2, a1, a0, F0, omega, phi0, x0, v0, t_span)
% SOLVE_LTI_TIME_INVARIANCE_SINUSOID Solves a2*x'' + a1*x' + a0*x = F0*cos(omega*t - phi0).
% Employs the Time-Invariance Principle and Control System Toolbox.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, lsim, bode, pole).
%
% % MBD Concept Mapping: Demonstrates Time-Invariance Principle: p(D)y(t-a) = q(t-a).
% %                      Validates that plant phase lag phi_plant depends strictly on omega,
% %                      and is invariant to the excitation initial phase phi0.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 17 Benchmark):
%   a2      - Coefficient of x'' (default: 1.0)
%   a1      - Coefficient of x' (default: 1.0)
%   a0      - Coefficient of x (default: 6.0)
%   F0      - Amplitude of sinusoidal forcing (default: 1.0)
%   omega   - Driving circular frequency in rad/s (default: 2.0)
%   phi0    - Input phase delay in radians (default: pi/3)
%   x0      - Initial position state x(0) (default: 0.0)
%   v0      - Initial velocity state x'(0) (default: 0.0)
%   t_span  - [t_start, t_end] simulation interval (default: [0, 16.0])
%
% Outputs:
%   sol_analytic - Struct containing function handles:
%                  .x_total, .x_steady_state, .x_transient, .u_input
%   sys_tf       - Transfer function object G(s) from Control System Toolbox
%   sim_data     - Struct containing lsim numerical verification results [t, y, u]

    arguments
        a2 (1,1) double {mustBePositive} = 1.0
        a1 (1,1) double {mustBeNonnegative} = 1.0
        a0 (1,1) double {mustBePositive} = 6.0
        F0 (1,1) double {mustBeNumeric} = 1.0
        omega (1,1) double {mustBePositive} = 2.0
        phi0 (1,1) double {mustBeNumeric} = pi / 3
        x0 (1,1) double {mustBeNumeric} = 0.0
        v0 (1,1) double {mustBeNumeric} = 0.0
        t_span (1,2) double {mustBeNumeric} = [0, 16.0]
    end

    %% 1. Plant Construction via Control System Toolbox
    sys_tf = tf(1, [a2, a1, a0]);

    %% 2. Steady-State Evaluation via Time Invariance & Complex Gain
    % Step A: Time delay parameter a = phi0 / omega
    delay_a = phi0 / omega;
    
    % Step B: Characteristic polynomial evaluation at s = 1i*omega
    p_val = a2 * (1i * omega)^2 + a1 * (1i * omega) + a0;
    
    if abs(p_val) < 1e-12
        error('Resonant condition: p(i*omega) vanishes.');
    end
    
    % Intrinsic Plant Gain and Phase Lag (Independent of phi0)
    gain_mag = 1 / abs(p_val);
    phi_plant_lag = angle(p_val); % Radians
    
    % Output steady-state amplitude
    A_ss = F0 * gain_mag;
    
    % Step C: Apply Time Invariance: x_p(t) = y_p(t - a)
    % y_p(t) = A_ss * cos(omega*t - phi_plant_lag)
    % x_p(t) = A_ss * cos(omega*(t - delay_a) - phi_plant_lag) = A_ss * cos(omega*t - phi0 - phi_plant_lag)
    total_phase_lag = phi0 + phi_plant_lag;
    x_p_handle = @(t) A_ss * cos(omega * t - total_phase_lag);
    v_p_handle = @(t) -A_ss * omega * sin(omega * t - total_phase_lag);

    x_p_0 = x_p_handle(0);
    v_p_0 = v_p_handle(0);

    %% 3. Transient Homogeneous Decomposition (Superposition II)
    disc = a1^2 - 4 * a2 * a0;
    x_h_0 = x0 - x_p_0;
    v_h_0 = v0 - v_p_0;

    if disc < -1e-12
        alpha = -a1 / (2 * a2);
        w_d = sqrt(4 * a2 * a0 - a1^2) / (2 * a2);
        c1 = x_h_0;
        c2 = (v_h_0 - alpha * c1) / w_d;
        x_h_handle = @(t) exp(alpha * t) .* (c1 * cos(w_d * t) + c2 * sin(w_d * t));
    elseif abs(disc) <= 1e-12
        r = -a1 / (2 * a2);
        c1 = x_h_0;
        c2 = v_h_0 - r * c1;
        x_h_handle = @(t) (c1 + c2 * t) .* exp(r * t);
    else
        r1 = (-a1 + sqrt(disc)) / (2 * a2);
        r2 = (-a1 - sqrt(disc)) / (2 * a2);
        c2 = (v_h_0 - r1 * x_h_0) / (r2 - r1);
        c1 = x_h_0 - c2;
        x_h_handle = @(t) c1 * exp(r1 * t) + c2 * exp(r2 * t);
    end

    x_total_handle = @(t) x_p_handle(t) + x_h_handle(t);

    sol_analytic.x_total = x_total_handle;
    sol_analytic.x_steady_state = x_p_handle;
    sol_analytic.x_transient = x_h_handle;
    sol_analytic.u_input = @(t) F0 * cos(omega * t - phi0);
    sol_analytic.gain = gain_mag;
    sol_analytic.phi_plant_lag_deg = phi_plant_lag * (180 / pi);
    sol_analytic.total_phase_lag_deg = total_phase_lag * (180 / pi);
    sol_analytic.delay_a = delay_a;

    %% 4. Numerical Verification via lsim
    t_eval = linspace(t_span(1), t_span(2), 1000);
    u_eval = sol_analytic.u_input(t_eval);
    
    if x0 == 0 && v0 == 0
        [y_sim, ~] = lsim(sys_tf, u_eval, t_eval);
    else
        ode_eq = @(t, X) [X(2); (F0*cos(omega*t - phi0) - a1*X(2) - a0*X(1))/a2];
        [~, X_ode] = ode45(ode_eq, t_eval, [x0; v0]);
        y_sim = X_ode(:, 1);
    end

    sim_data.t = t_eval;
    sim_data.u = u_eval;
    sim_data.y = y_sim;

    max_error = max(abs(x_total_handle(t_eval)' - y_sim(:)));
    sol_analytic.verification_max_error = max_error;

    %% 5. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Time-Invariance LTI Solver: MIT 18.03 Lec 17     \n');
    fprintf('====================================================\n');
    fprintf('ODE: %.2f*x'''' + %.2f*x'' + %.2f*x = %.2f*cos(%.2f*t - %.4f)\n', ...
            a2, a1, a0, F0, omega, phi0);
    fprintf('Input Phase Shift (phi0): %.4f rad (%.2f deg) -> Delay a = %.4f s\n', ...
            phi0, phi0*180/pi, delay_a);
    fprintf('Plant Operator p(i*w)   : %+.4f %+.4fi\n', real(p_val), imag(p_val));
    fprintf('Intrinsic Plant Gain    : %.4f (1 / |p(i*w)|)\n', gain_mag);
    fprintf('Intrinsic Plant Lag     : %.4f rad (%.2f deg)\n', ...
            phi_plant_lag, sol_analytic.phi_plant_lag_deg);
    fprintf('Steady Output Amplitude : %.4f\n', A_ss);
    fprintf('Total Output Lag        : %.4f rad (%.2f deg)\n', ...
            total_phase_lag, sol_analytic.total_phase_lag_deg);
    fprintf('Verification Max Error  : %e\n', max_error);
    fprintf('====================================================\n\n');

    %% 6. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 950, 500]);
        hold on; grid on; box on;

        plot(t_eval, u_eval, 'k:', 'LineWidth', 1.2, ...
             'DisplayName', sprintf('Input Signal u(t) (\\phi_0 = %.1f^\\circ)', phi0*180/pi));
        plot(t_eval, x_total_handle(t_eval), 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Total Analytical Response x(t)');
        plot(t_eval, x_p_handle(t_eval), 'r--', 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Steady x_p(t) (Amp = %.3f, Lag = %.1f^\\circ)', ...
                                    A_ss, sol_analytic.total_phase_lag_deg));
        plot(t_eval(1:40:end), y_sim(1:40:end), 'mo', 'MarkerSize', 5, ...
             'DisplayName', 'lsim Numerical Check');

        xlabel('Time t (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('State Position x(t)', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('Time Invariance: Delay a = %.3fs, Plant Lag = %.1f^\\circ, Total Lag = %.1f^\\circ', ...
              delay_a, sol_analytic.phi_plant_lag_deg, sol_analytic.total_phase_lag_deg), ...
              'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end