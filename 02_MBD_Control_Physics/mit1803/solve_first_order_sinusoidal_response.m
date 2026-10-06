function [sol_analytic, sys_metrics] = solve_first_order_sinusoidal_response(k, F0, omega, phi_in, x0, t_span)
% SOLVE_FIRST_ORDER_SINUSOIDAL_RESPONSE Solves x' + k*x = F0*cos(omega*t - phi_in).
%
% % Required Toolboxes: 
% %   - Base MATLAB (Zero external dependency).
%
% % MBD Concept Mapping: Demonstrates frequency response G(j*omega) = 1/(k + j*omega).
% %                      Computes sinusoidal steady-state tracking gain and phase lag,
% %                      and separates decaying transient dynamics.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 7 Benchmark):
%   k        - Plant pole parameter (x' + k*x), default: 2.0
%   F0       - Amplitude of sinusoidal excitation, default: 2.0
%   omega    - Circular excitation frequency (rad/s), default: 1.0
%   phi_in   - Input phase lag in radians, default: 0.0
%   x0       - Initial state condition x(0), default: 0.0
%   t_span   - [t_start, t_end] simulation interval, default: [0, 10.0]
%
% Outputs:
%   sol_analytic - Struct containing closed-form function handles:
%                  .x_total, .x_steady_state, .x_transient
%   sys_metrics  - Struct containing gains, phase delays, time constant, and error

    arguments
        k (1,1) double {mustBePositive} = 2.0
        F0 (1,1) double {mustBeNumeric} = 2.0
        omega (1,1) double {mustBePositive} = 1.0
        phi_in (1,1) double {mustBeNumeric} = 0.0
        x0 (1,1) double {mustBeNumeric} = 0.0
        t_span (1,2) double {mustBeNumeric} = [0, 10.0]
    end

    %% 1. Analytical Evaluation via Complex Replacement (Phasor Method)
    % Frequency response of G(s) = 1 / (s + k) at s = j*omega
    G_jw = 1 / (k + 1i * omega);
    
    % Input phasor in polar representation: U_in = F0 * exp(-1i * phi_in)
    U_in = F0 * exp(-1i * phi_in);
    
    % Output complex phasor: Z_out = G(j*omega) * U_in
    Z_out = G_jw * U_in;
    
    A_ss = abs(Z_out);                    % Output sinusoidal amplitude
    phi_total_lag = -angle(Z_out);         % Output phase lag relative to cos(omega*t)
    
    % Steady-state particular solution: x_p(t) = Re(Z_out * exp(1i * omega * t))
    x_p_handle = @(t) real(Z_out * exp(1i * omega * t));
    
    % Transient homogeneous decay: x_h(t) = C * exp(-k*t)
    % Enforcing initial state condition x(0) = x_p(0) + C -> C = x0 - x_p(0)
    C_const = x0 - real(Z_out);
    x_h_handle = @(t) C_const * exp(-k * t);
    
    % Full closed-form analytical solution handle
    x_total_handle = @(t) x_p_handle(t) + x_h_handle(t);

    sol_analytic.x_total = x_total_handle;
    sol_analytic.x_steady_state = x_p_handle;
    sol_analytic.x_transient = x_h_handle;
    sol_analytic.C_constant = C_const;

    %% 2. Performance Metrics Packaging
    sys_metrics.time_constant = 1 / k;
    sys_metrics.system_gain = abs(G_jw);
    sys_metrics.gain_dB = 20 * log10(sys_metrics.system_gain);
    sys_metrics.phase_lag_rad = phi_total_lag;
    sys_metrics.phase_lag_deg = phi_total_lag * (180 / pi);
    sys_metrics.steady_amplitude = A_ss;

    %% 3. Numerical Cross-Check via ode45
    t_eval = linspace(t_span(1), t_span(2), 600);
    u_fun = @(t) F0 * cos(omega * t - phi_in);
    ode_fun = @(t, x) -k * x + u_fun(t);
    
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    [t_num, x_num] = ode45(ode_fun, t_eval, x0, opts);
    
    x_ana_eval = x_total_handle(t_num);
    max_error = max(abs(x_ana_eval - x_num));
    sys_metrics.verification_max_error = max_error;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   First-Order Sinusoidal Solver: MIT 18.03 Lec 7   \n');
    fprintf('====================================================\n');
    fprintf('Differential Equation   : dx/dt + %.2f*x = %.2f*cos(%.2f*t - %.2f)\n', ...
            k, F0, omega, phi_in);
    fprintf('System Pole             : s = -%.2f (Time Constant tau = %.4f s)\n', k, sys_metrics.time_constant);
    fprintf('Frequency Response G(jw): %+.4f %+.4fi\n', real(G_jw), imag(G_jw));
    fprintf('Steady Gain Magnitude   : %.4f (%.2f dB)\n', sys_metrics.system_gain, sys_metrics.gain_dB);
    fprintf('Steady-State Amplitude  : %.4f\n', sys_metrics.steady_amplitude);
    fprintf('Steady-State Phase Lag  : %.4f rad (%.2f degrees)\n', ...
            sys_metrics.phase_lag_rad, sys_metrics.phase_lag_deg);
    fprintf('Transient Constant C    : %.4f (Decay mode: e^{-%.2f t})\n', C_const, k);
    fprintf('Verification Max Error  : %e\n', max_error);
    fprintf('====================================================\n\n');

    %% 5. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 900, 520]);
        hold on; grid on; box on;

        plot(t_eval, u_fun(t_eval), 'k:', 'LineWidth', 1.2, ...
             'DisplayName', sprintf('Input Signal u(t) (Amp = %.1f)', F0));
        plot(t_eval, x_total_handle(t_eval), 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Analytic Total Response x(t)');
        plot(t_eval, x_p_handle(t_eval), 'r--', 'LineWidth', 1.6, ...
             'DisplayName', sprintf('Steady-State x_p(t) (Amp = %.3f)', A_ss));
        plot(t_eval, x_h_handle(t_eval), 'g-.', 'LineWidth', 1.4, ...
             'DisplayName', 'Decaying Transient x_h(t)');
        plot(t_num(1:25:end), x_num(1:25:end), 'mo', 'MarkerSize', 5, ...
             'DisplayName', 'ode45 Numerical Check');

        xlabel('Time t (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('State Variable x(t)', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('LTI Sinusoidal Response: Gain = %.3f, Phase Lag = %.1f^\\circ', ...
              sys_metrics.system_gain, sys_metrics.phase_lag_deg), ...
              'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'southeast', 'FontSize', 9);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end