function [ident_results, sys_tf] = estimate_damping_ratio_half_amplitude(T_10_cycles, n_50_cycles, x0, num_periods_plot)
% ESTIMATE_DAMPING_RATIO_HALF_AMPLITUDE Identifies 2nd-order parameters via the 50% decay rule.
% Governing Rule (Kim Vandiver, MIT 18.03 Lec 18): zeta ~ 0.11 / n_50%.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, damp, initial, ss, pole).
%
% % MBD Concept Mapping: Encapsulates standard 2nd-order plant: G(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2).
% %                      Implements modal parameter identification via logarithmic decrement,
% %                      extracting continuous damping ratio zeta and natural frequency wn.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 18 Demo):
%   T_10_cycles       - Measured elapsed time for 10 complete cycles in seconds (default: 14.2)
%   n_50_cycles       - Number of cycles to decay to 50% amplitude (default: 4.0)
%   x0                - Initial release displacement (default: 1.0)
%   num_periods_plot  - Number of periods to simulate for validation plot (default: 8.0)
%
% Outputs:
%   ident_results - Struct containing identified parameters:
%                   .zeta_approx, .zeta_exact, .omega_d, .omega_n, .period, .poles
%   sys_tf        - Standard second-order transfer function object G(s)

    arguments
        T_10_cycles (1,1) double {mustBePositive} = 14.2
        n_50_cycles (1,1) double {mustBePositive} = 4.0
        x0 (1,1) double {mustBeNumeric} = 1.0
        num_periods_plot (1,1) double {mustBePositive} = 8.0
    end

    %% 1. Frequency and Period Evaluation
    P_val = T_10_cycles / 10.0;     % Period per cycle (seconds)
    f_d = 1.0 / P_val;              % Damped frequency in Hz
    omega_d = 2.0 * pi * f_d;       % Damped circular frequency (rad/s)

    %% 2. Damping Ratio Identification
    % Vandiver Quick Engineering Formula: zeta ~ 0.11 / n_50
    zeta_approx = 0.11 / n_50_cycles;

    % Exact Logarithmic Decrement Derivation:
    % exp(-n_50 * zeta * wn * P) = 0.5 -> delta_50 = ln(2) / n_50
    delta_50 = log(2.0) / n_50_cycles;
    zeta_exact = delta_50 / sqrt((2.0 * pi)^2 + delta_50^2);

    % Undamped Natural Frequency
    omega_n = omega_d / sqrt(1.0 - zeta_exact^2);

    %% 3. Control System Toolbox LTI Model Construction
    num_tf = omega_n^2;
    den_tf = [1.0, 2.0 * zeta_exact * omega_n, omega_n^2];
    sys_tf = tf(num_tf, den_tf);
    sys_poles = pole(sys_tf);

    ident_results.period = P_val;
    ident_results.f_d_hz = f_d;
    ident_results.omega_d = omega_d;
    ident_results.omega_n = omega_n;
    ident_results.zeta_approx = zeta_approx;
    ident_results.zeta_exact = zeta_exact;
    ident_results.damping_percent = zeta_exact * 100.0;
    ident_results.poles = sys_poles;

    %% 4. Time-Domain Response Simulation via Control System Toolbox
    t_end = num_periods_plot * P_val;
    t_eval = linspace(0, t_end, 1000);
    
    % Realize state-space for zero-input initial displacement simulation
    A_mat = [0, 1; -omega_n^2, -2*zeta_exact*omega_n];
    B_mat = [0; 1];
    C_mat = [1, 0];
    D_mat = 0;
    sys_ss = ss(A_mat, B_mat, C_mat, D_mat);
    
    [y_sim, t_out] = initial(sys_ss, [x0; 0], t_eval);

    %% 5. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Vibration Damping Estimator: MIT 18.03 Lec 18    \n');
    fprintf('====================================================\n');
    fprintf('Elapsed Time for 10 Cycles: %.4f s -> Period P = %.4f s\n', T_10_cycles, P_val);
    fprintf('Damped Frequency (w_d)    : %.4f rad/s (%.3f Hz)\n', omega_d, f_d);
    fprintf('Cycles to 50%% Decay (n_50): %.2f cycles\n', n_50_cycles);
    fprintf('Vandiver Approx (zeta)    : %.4f (%.2f%% critical)\n', zeta_approx, zeta_approx * 100);
    fprintf('Exact Identified (zeta)   : %.4f (%.2f%% critical)\n', zeta_exact, zeta_exact * 100);
    fprintf('Natural Frequency (w_n)   : %.4f rad/s\n', omega_n);
    fprintf('\n--- Control System Toolbox Modal Poles ---\n');
    damp(sys_tf);
    fprintf('====================================================\n\n');

    %% 6. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 900, 480]);
        hold on; grid on; box on;

        % Plot simulated trajectory
        plot(t_out, y_sim, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Identified Response x(t)');

        % Plot decay envelopes
        env_up = x0 * exp(-zeta_exact * omega_n * t_eval);
        plot(t_eval, env_up, 'r--', 'LineWidth', 1.4, 'DisplayName', 'Decay Envelope \pm e^{-\zeta \omega_n t}');
        plot(t_eval, -env_up, 'r--', 'LineWidth', 1.4, 'HandleVisibility', 'off');

        % Highlight 50% decay point
        t_50 = n_50_cycles * P_val;
        yline(0.5 * x0, 'k:', 'LineWidth', 1.2, 'DisplayName', '50% Amplitude Level');
        plot(t_50, 0.5 * x0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7, ...
             'DisplayName', sprintf('Half-Amplitude Point (t = %.2fs, n = %.1f)', t_50, n_50_cycles));

        xlabel('Time t (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('State Displacement x(t)', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('Vandiver Damping Rule Identification: $\\zeta = %.2f\\%%$ ($n_{50\\%%} = %.1f$)', ...
              ident_results.damping_percent, n_50_cycles), 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter','latex');
        legend('Location', 'northeast', 'FontSize', 9);
        xlim([0, t_end]);
        hold off;
    end
end