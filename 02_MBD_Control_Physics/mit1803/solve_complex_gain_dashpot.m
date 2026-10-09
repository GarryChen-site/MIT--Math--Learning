function [sol_analytic, sys_metrics] = solve_complex_gain_dashpot(m, b, k, A, omega, x0, v0, t_span)
% SOLVE_COMPLEX_GAIN_DASHPOT Solves dashpot-driven ODE: m*x'' + b*x' + k*x = b*y'
% where y(t) = A*cos(omega*t) is the base displacement excitation.
%
% % Required Toolboxes: 
% %   - Base MATLAB (Zero external dependency).
%
% % MBD Concept Mapping: Evaluates continuous Complex Gain H(w) = b*i*w / p(i*w).
% %                      Demonstrates band-pass transfer operator G(s) = b*s / (m*s^2 + b*s + k)
% %                      governing quarter-car suspension and base vibration isolation.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 14 Benchmark):
%   m      - Mass coefficient (default: 1.0)
%   b      - Damping coefficient (default: 1.0)
%   k      - Spring stiffness coefficient (default: 2.0)
%   A      - Input base displacement amplitude (default: 1.0)
%   omega  - Input driving circular frequency (default: 1.0)
%   x0     - Initial position state x(0) (default: 0.0)
%   v0     - Initial velocity state x'(0) (default: 0.0)
%   t_span - [t_start, t_end] simulation interval (default: [0, 16.0])
%
% Outputs:
%   sol_analytic - Struct containing function handles:
%                  .x_total, .x_steady_state, .x_transient, .y_input
%   sys_metrics  - Struct containing Complex Gain H, magnitude, phase lag, and errors

    arguments
        m (1,1) double {mustBePositive} = 1.0
        b (1,1) double {mustBeNonnegative} = 1.0
        k (1,1) double {mustBePositive} = 2.0
        A (1,1) double {mustBeNumeric} = 1.0
        omega (1,1) double {mustBePositive} = 1.0
        x0 (1,1) double {mustBeNumeric} = 0.0
        v0 (1,1) double {mustBeNumeric} = 0.0
        t_span (1,2) double {mustBeNumeric} = [0, 16.0]
    end

    %% 1. Analytical Complex Gain Derivation
    % Characteristic polynomial evaluated at s = 1i*omega:
    p_val = m * (1i * omega)^2 + b * (1i * omega) + k;
    
    if abs(p_val) < 1e-12
        error('Resonance condition: p(i*omega) vanishes.');
    end
    
    % Complex Gain for dashpot driving: H(w) = b*i*w / p(i*w)
    H_complex = (1i * b * omega) / p_val;
    
    gain_mag = abs(H_complex);
    phi_lag = -angle(H_complex); % Phase delay relative to cos(omega*t)
    
    % Steady-state particular solution: x_p(t) = A * |H| * cos(omega*t - phi)
    x_p_handle = @(t) A * gain_mag * cos(omega * t - phi_lag);
    v_p_handle = @(t) -A * gain_mag * omega * sin(omega * t - phi_lag);
    
    x_p_0 = x_p_handle(0);
    v_p_0 = v_p_handle(0);

    %% 2. Transient Homogeneous Decomposition
    disc = b^2 - 4 * m * k;
    x_h_0 = x0 - x_p_0;
    v_h_0 = v0 - v_p_0;
    
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

    sol_analytic.y_input = @(t) A * cos(omega * t);
    sol_analytic.x_steady_state = x_p_handle;
    sol_analytic.x_transient = x_h_handle;
    sol_analytic.x_total = @(t) x_p_handle(t) + x_h_handle(t);

    %% 3. Metrics Packaging
    sys_metrics.H_complex = H_complex;
    sys_metrics.gain = gain_mag;
    sys_metrics.phase_lag_rad = phi_lag;
    sys_metrics.phase_lag_deg = phi_lag * (180 / pi);
    sys_metrics.p_eval = p_val;

    %% 4. Numerical Verification via ode45
    % Forcing term: b*y'(t) = -b*A*omega*sin(omega*t)
    ode_fun = @(t, X) [X(2); (-b*A*omega*sin(omega*t) - b*X(2) - k*X(1)) / m];
    t_eval = linspace(t_span(1), t_span(2), 600);
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    
    [t_num, X_num] = ode45(ode_fun, t_eval, [x0; v0], opts);
    x_analytic_eval = sol_analytic.x_total(t_num);
    max_error = max(abs(x_analytic_eval - X_num(:, 1)));
    sys_metrics.verification_max_error = max_error;

    %% 5. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Dashpot-Driven Solver: MIT 18.03 Lecture 14      \n');
    fprintf('====================================================\n');
    fprintf('Physical System         : %.2f*x'''' + %.2f*x'' + %.2f*x = %.2f*y''\n', m, b, k, b);
    fprintf('Input Excitation y(t)   : %.2f*cos(%.2f*t)\n', A, omega);
    fprintf('Characteristic p(iw)    : %+.4f %+.4fi\n', real(p_val), imag(p_val));
    fprintf('Complex Gain H(w)       : %+.4f %+.4fi\n', real(H_complex), imag(H_complex));
    fprintf('Gain Magnitude |H(w)|   : %.4f (Analytical: sqrt(2)/2 = 0.7071)\n', gain_mag);
    fprintf('Phase Lag phi = -Arg(H) : %+.4f rad (%+.1f degrees)\n', phi_lag, phi_lag * 180/pi);
    fprintf('Verification Max Error  : %e\n', max_error);
    fprintf('====================================================\n\n');

    %% 6. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 950, 500]);
        hold on; grid on; box on;

        plot(t_eval, sol_analytic.y_input(t_eval), 'k:', 'LineWidth', 1.2, ...
             'DisplayName', sprintf('Base Input y(t) = %.1fcos(%.1ft)', A, omega));
        plot(t_eval, sol_analytic.x_total(t_eval), 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Total Mass Position x(t)');
        plot(t_eval, sol_analytic.x_steady_state(t_eval), 'r--', 'LineWidth', 1.6, ...
             'DisplayName', sprintf('Steady Output x_p(t) (Gain=%.3f, Lag=%+.1f^\\circ)', ...
                                    gain_mag, sys_metrics.phase_lag_deg));
        plot(t_eval, sol_analytic.x_transient(t_eval), 'g-.', 'LineWidth', 1.4, ...
             'DisplayName', 'Transient Decay x_h(t)');
        plot(t_num(1:25:end), X_num(1:25:end, 1), 'mo', 'MarkerSize', 5, ...
             'DisplayName', 'ode45 Numerical Check');

        xlabel('Time t (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Displacement', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('Dashpot Driving Response: |H| = %.3f, Phase Lag \\phi = %+.1f^\\circ', ...
              gain_mag, sys_metrics.phase_lag_deg), 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end