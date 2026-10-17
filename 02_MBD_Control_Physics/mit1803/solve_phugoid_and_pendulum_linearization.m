function [aero_metrics, pendulum_metrics, sys_models] = solve_phugoid_and_pendulum_linearization(v0_speed, thrust_weight_ratio, L_pend, b_pend, pend_ic)
% SOLVE_PHUGOID_AND_PENDULUM_LINEARIZATION Analyzes Airplane Phugoid Dynamics & Nonlinear Pendulum.
% Computes Phugoid natural frequency wn = sqrt(2)*g/v0 and damping ratio zeta = (F/mg)/sqrt(2),
% and simulates the nonlinear pendulum x' = y, y' = -(g/L)*sin(x) - b*y across multiple basins.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (ss, initial, pole, damp).
%
% % MBD Concept Mapping: Connects linearized state-space flight dynamics to standard 2nd-order
% %                      specifications (wn, zeta), and tracks nonlinear basin-of-attraction capture.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 38 Boeing 747 & Pendulum Benchmarks):
%   v0_speed            - Equilibrium cruise speed v0 in m/s (default: 260.0 for Boeing 747)
%   thrust_weight_ratio - Dimensionless Thrust-to-Weight ratio F/(m*g) (default: 0.27 for Boeing 747)
%   L_pend              - Length of nonlinear pendulum rod in meters (default: 2.4525 -> k = g/L = 4.0)
%   b_pend              - Pendulum damping coefficient c/m in 1/s (default: 0.4)
%   pend_ic             - Initial pendulum state [theta(0); theta'(0)] (default: [-pi; 5.0])
%
% Outputs:
%   aero_metrics     - Struct with Phugoid matrix A, wn, zeta, period_undamped, and period_damped
%   pendulum_metrics - Struct with k, Downward/Upward Jacobians, final attractor 2*n*pi, and trajectory
%   sys_models       - Struct containing Control System Toolbox 'ss' objects for Phugoid and Pendulum

    arguments
        v0_speed (1,1) double {mustBePositive} = 260.0
        thrust_weight_ratio (1,1) double {mustBeNonnegative} = 0.27
        L_pend (1,1) double {mustBePositive} = 2.4525
        b_pend (1,1) double {mustBeNonnegative} = 0.4
        pend_ic (2,1) double {mustBeNumeric} = [-pi; 5.0]
    end

    g_acc = 9.81;

    %% 1. Airplane Phugoid State-Space Model & Modal Analysis (Section [2])
    % State equation: d/dt [u; w] = [ -2*F/(m*v0), -g/v0 ; 2*g/v0, 0 ] * [u; w]
    A_phugoid = [-2.0 * thrust_weight_ratio * (g_acc / v0_speed), -g_acc / v0_speed;
                  2.0 * g_acc / v0_speed,                          0.0             ];

    sys_phugoid = ss(A_phugoid, zeros(2, 1), eye(2), zeros(2, 1));
    [wn_vec, zeta_vec, poles_phugoid] = damp(sys_phugoid);

    wn_phugoid   = wn_vec(1);
    zeta_phugoid = zeta_vec(1);
    T_undamped   = 2.0 * pi / wn_phugoid; % Lecture formula: ~ 0.45 * v0
    if zeta_phugoid < 1.0
        T_damped = 2.0 * pi / (wn_phugoid * sqrt(1.0 - zeta_phugoid^2));
    else
        T_damped = Inf;
    end

    aero_metrics.A_phugoid = A_phugoid;
    aero_metrics.v0 = v0_speed;
    aero_metrics.thrust_to_weight = thrust_weight_ratio;
    aero_metrics.omega_n = wn_phugoid;
    aero_metrics.zeta = zeta_phugoid;
    aero_metrics.poles = poles_phugoid;
    aero_metrics.period_undamped_s = T_undamped;
    aero_metrics.period_damped_s = T_damped;

    %% 2. Nonlinear Pendulum Linearization & Multi-Basin Simulation (Section [1])
    k_pend = g_acc / L_pend;
    J_attractor = [0, 1; -k_pend, -b_pend]; % At theta = 2*n*pi
    J_saddle    = [0, 1; +k_pend, -b_pend]; % At theta = (2n+1)*pi

    sys_pend_attractor = ss(J_attractor, zeros(2, 1), eye(2), zeros(2, 1));
    sys_pend_saddle    = ss(J_saddle,    zeros(2, 1), eye(2), zeros(2, 1));

    % Integrate full nonlinear pendulum until it settles into a basin attractor
    t_pend = linspace(0, 20.0, 1200);
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    ode_pend = @(t, U) [U(2); -k_pend * sin(U(1)) - b_pend * U(2)];
    [~, U_pend_traj] = ode45(ode_pend, t_pend, pend_ic, opts);

    % Identify which 2*n*pi attractor captured the pendulum
    final_theta = U_pend_traj(end, 1);
    basin_index_n = round(final_theta / (2 * pi));
    captured_attractor_rad = basin_index_n * 2 * pi;

    pendulum_metrics.k = k_pend;
    pendulum_metrics.b = b_pend;
    pendulum_metrics.J_attractor = J_attractor;
    pendulum_metrics.J_saddle = J_saddle;
    pendulum_metrics.attractor_poles = pole(sys_pend_attractor);
    pendulum_metrics.saddle_poles = pole(sys_pend_saddle);
    pendulum_metrics.captured_basin_multiple_of_2pi = basin_index_n;
    pendulum_metrics.captured_attractor_rad = captured_attractor_rad;
    pendulum_metrics.t = t_pend;
    pendulum_metrics.U_traj = U_pend_traj;

    sys_models.sys_phugoid = sys_phugoid;
    sys_models.sys_pend_attractor = sys_pend_attractor;
    sys_models.sys_pend_saddle = sys_pend_saddle;

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Phugoid & Pendulum Solver: MIT 18.03 Lecture 38  \n');
    fprintf('====================================================\n');
    fprintf('--- [1] Airplane Phugoid Oscillation ---\n');
    fprintf('Cruise Speed (v0)       : %.2f m/s\n', v0_speed);
    fprintf('Thrust-to-Weight (F/mg) : %.4f (Dimensionless)\n', thrust_weight_ratio);
    fprintf('Natural Frequency (w_n) : %.5f rad/s (sqrt(2)*g/v0)\n', wn_phugoid);
    fprintf('Damping Ratio (zeta)    : %.4f ((F/mg)/sqrt(2))\n', zeta_phugoid);
    fprintf('Phugoid Period (2pi/wn) : %.2f seconds (~0.45 * v0)\n', T_undamped);
    fprintf('\n--- [2] Nonlinear Damped Pendulum ---\n');
    fprintf('Parameters              : k = g/L = %.2f s^-2, b = %.2f s^-1\n', k_pend, b_pend);
    fprintf('Initial State           : theta(0) = %.2f*pi rad, theta''(0) = %.2f rad/s\n', ...
            pend_ic(1)/pi, pend_ic(2));
    fprintf('Captured Attractor Basin: theta_eq = %d * (2*pi) = %.4f rad\n', ...
            basin_index_n, captured_attractor_rad);
    fprintf('====================================================\n\n');

    %% 4. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [100, 100, 1140, 460]);

        % Left Subplot: Phugoid Velocity Perturbation Phase Plane & Time Trace
        subplot(1, 2, 1);
        hold on; grid on; box on;
        t_ph = linspace(0, 3.5 * T_undamped, 800);
        [Y_ph, ~] = initial(sys_phugoid, [0; 10.0], t_ph); % 10 m/s vertical downdraft
        plot(t_ph, Y_ph(:, 1), 'b-', 'LineWidth', 2.0, 'DisplayName', 'Horizontal Perturbation u(t)');
        plot(t_ph, Y_ph(:, 2), 'r--', 'LineWidth', 1.8, 'DisplayName', 'Vertical Perturbation w(t)');
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Velocity Perturbation (m/s)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Phugoid Mode (v_0=%.0f m/s, T=%.1fs, \\zeta=%.2f)', ...
              v0_speed, T_undamped, zeta_phugoid), 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([0, t_ph(end)]);

        % Right Subplot: Pendulum Over-the-Top Trajectory & Basin Capture
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(U_pend_traj(:, 1)/pi, U_pend_traj(:, 2), 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Pendulum Trajectory (\theta/\pi, \theta'')');
        plot(pend_ic(1)/pi, pend_ic(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, ...
             'DisplayName', 'Initial State');
        plot(2*basin_index_n, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 11, ...
             'DisplayName', sprintf('Captured Attractor (%d\\pi, 0)', 2*basin_index_n));

        % Mark intermediate saddles passed along the way
        saddles_n = -3:2:5;
        plot(saddles_n, zeros(size(saddles_n)), 'kd', 'MarkerFaceColor', 'y', ...
             'MarkerSize', 6, 'DisplayName', 'Upward Saddles ((2k+1)\pi)');

        xlabel('Angle \theta / \pi (multiples of \pi rad)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Angular Velocity \theta'' (rad/s)', 'FontSize', 10, 'FontWeight', 'bold');
        title('Nonlinear Pendulum: Full Rotations & Basin Capture', 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        hold off;
    end
end