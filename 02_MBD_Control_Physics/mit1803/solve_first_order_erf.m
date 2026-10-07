function [sol_analytic, sys_metrics] = solve_first_order_erf(k, B, omega, theta_in, x0, t_span)
% SOLVE_FIRST_ORDER_ERF Solves dx/dt + k*x = B*cos(omega*t + theta_in) via ERF.
%
% % Required Toolboxes: 
% %   - Base MATLAB (Zero external dependency).
%
% % MBD Concept Mapping: Demonstrates frequency-domain evaluation of LTI transfer
% %                      operator G(s) = 1/P(s) via the Exponential Response Formula (ERF).
% %                      Computes sinusoidal steady-state amplitude scaling and phase shift.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 9 Benchmark):
%   k        - ODE state coefficient (dx/dt + k*x), default: -1.0 (Unstable pole s = 1)
%   B        - Amplitude of sinusoidal input, default: 3.0
%   omega    - Input circular frequency (rad/s), default: 1.0
%   theta_in - Input phase angle in radians, default: 0.0
%   x0       - Initial state condition x(0), default: -1.5 (Matches periodic orbit)
%   t_span   - [t_start, t_end] simulation interval, default: [0, 6.0]
%
% Outputs:
%   sol_analytic - Struct containing function handles:
%                  .x_total, .x_steady_state, .x_homogeneous
%   sys_metrics  - Struct containing ERF evaluation, amplitude, phase, and validation errors

    arguments
        k (1,1) double {mustBeNumeric} = -1.0
        B (1,1) double {mustBeNumeric} = 3.0
        omega (1,1) double {mustBePositive} = 1.0
        theta_in (1,1) double {mustBeNumeric} = 0.0
        x0 (1,1) double {mustBeNumeric} = -1.5
        t_span (1,2) double {mustBeNumeric} = [0, 6.0]
    end

    %% 1. Analytical Solution via Exponential Response Formula (ERF)
    % Characteristic polynomial: P(s) = s + k
    % Frequency response evaluated at input frequency: P(j*omega) = 1i*omega + k
    P_val = 1i * omega + k;
    
    if abs(P_val) < 1e-12
        error('Resonance condition: P(j*omega) = 0. ERF denominator vanishes.');
    end
    
    % Input phasor: U = B * exp(1i * theta_in)
    U_in = B * exp(1i * theta_in);
    
    % ERF complex amplitude: Z_p = U / P(j*omega)
    Z_p = U_in / P_val;
    
    A_ss = abs(Z_p);               % Output steady-state amplitude
    phi_ss = angle(Z_p);           % Output phase angle
    
    % Steady-state particular solution: x_p(t) = Re(Z_p * exp(j*omega*t))
    x_p_handle = @(t) real(Z_p * exp(1i * omega * t));
    
    % Homogeneous component: x_h(t) = c * exp(-k*t)
    % Initial condition match: x(0) = x_p(0) + c -> c = x0 - x_p(0)
    c_const = x0 - real(Z_p);
    x_h_handle = @(t) c_const * exp(-k * t);
    
    % Full analytical solution handle
    x_total_handle = @(t) x_p_handle(t) + x_h_handle(t);

    sol_analytic.x_total = x_total_handle;
    sol_analytic.x_steady_state = x_p_handle;
    sol_analytic.x_homogeneous = x_h_handle;
    sol_analytic.c_constant = c_const;

    %% 2. Performance Metrics Packaging
    sys_metrics.pole = -k;
    sys_metrics.P_complex = P_val;
    sys_metrics.Z_complex = Z_p;
    sys_metrics.steady_amplitude = A_ss;
    sys_metrics.steady_phase_rad = phi_ss;
    sys_metrics.steady_phase_deg = phi_ss * (180 / pi);

    %% 3. Numerical Cross-Check via ode45
    t_eval = linspace(t_span(1), t_span(2), 600);
    u_fun = @(t) B * cos(omega * t + theta_in);
    ode_fun = @(t, x) -k * x + u_fun(t);
    
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [t_num, x_num] = ode45(ode_fun, t_eval, x0, opts);
    
    x_ana_eval = x_total_handle(t_num);
    max_error = max(abs(x_ana_eval - x_num));
    sys_metrics.verification_max_error = max_error;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   First-Order ERF Solver: MIT 18.03 Lecture 9      \n');
    fprintf('====================================================\n');
    fprintf('Differential Equation   : dx/dt %+.2f*x = %.2f*cos(%.2f*t %+.2f)\n', ...
            k, B, omega, theta_in);
    fprintf('System Pole (s)         : %+.4f\n', sys_metrics.pole);
    fprintf('Polynomial P(j*omega)   : %+.4f %+.4fi\n', real(P_val), imag(P_val));
    fprintf('ERF Phasor Z_p          : %+.4f %+.4fi\n', real(Z_p), imag(Z_p));
    fprintf('Steady Amplitude (A)    : %.4f\n', sys_metrics.steady_amplitude);
    fprintf('Steady Phase Angle      : %.4f rad (%.2f deg)\n', ...
            sys_metrics.steady_phase_rad, sys_metrics.steady_phase_deg);
    fprintf('Homogeneous Constant (c): %.4f\n', c_const);
    fprintf('Verification Max Error  : %e\n', max_error);
    fprintf('====================================================\n\n');

    %% 5. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 900, 500]);
        hold on; grid on; box on;

        plot(t_eval, u_fun(t_eval), 'k:', 'LineWidth', 1.2, ...
             'DisplayName', sprintf('Input u(t) = %.1fcos(%.1ft)', B, omega));
        plot(t_eval, x_total_handle(t_eval), 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Analytic Total Solution x(t)');
        plot(t_eval, x_p_handle(t_eval), 'r--', 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Steady-State ERF x_p(t) (Amp=%.3f)', A_ss));
        plot(t_num(1:25:end), x_num(1:25:end), 'mo', 'MarkerSize', 5, ...
             'DisplayName', 'ode45 Numerical Check');

        xlabel('Time t (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('State Variable x(t)', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('ERF Response: $dx/dt %+.1fx = %.1fcos(%.1ft), x(0)=%.2f$', ...
              k, B, omega, x0), 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter','latex');
        legend('Location', 'best', 'FontSize', 9);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end