function [vdp_report, lorenz_report, sys_models] = solve_nonlinear_limit_cycle_and_chaos(a_vdp, r_lorenz, t_span_vdp, t_span_lor)
% SOLVE_NONLINEAR_LIMIT_CYCLE_AND_CHAOS Analyzes van der Pol Limit Cycles & 3D Lorenz Chaos.
% Extracts local linearized 'ss' pole spectra at equilibria and identifies steady-state
% limit cycle amplitude and period for the van der Pol oscillator.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (ss, pole, damp).
%
% % MBD Concept Mapping: Contrasts local Jacobian pole analysis against global nonlinear attractors
% %                      (2D Limit Cycles and 3D Chaotic Strange Attractors).
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 39 Benchmarks):
%   a_vdp      - Nonlinear damping parameter a > 0 in x'' + a*(x^2 - 1)*x' + x = 0 (default: 1.5)
%   r_lorenz   - Lorenz Rayleigh parameter r in y' = -x*z + r*x - y (default: 28.0)
%   t_span_vdp - [t_start, t_end] for van der Pol limit cycle extraction (default: [0, 40.0])
%   t_span_lor - [t_start, t_end] for 3D Lorenz chaotic simulation (default: [0, 30.0])
%
% Outputs:
%   vdp_report    - Struct with origin Jacobian, poles, limit cycle amplitude, and period
%   lorenz_report - Struct with equilibria coordinates, 3x3 Jacobians, and 3D pole spectra
%   sys_models    - Struct containing Control System Toolbox 'ss' objects at equilibria

    arguments
        a_vdp (1,1) double {mustBePositive} = 1.5
        r_lorenz (1,1) double {mustBePositive} = 28.0
        t_span_vdp (1,2) double {mustBeNumeric} = [0.0, 40.0]
        t_span_lor (1,2) double {mustBeNumeric} = [0.0, 30.0]
    end

    %% 1. van der Pol Oscillator Analysis & Limit Cycle Identification (Section [3])
    % Companion system: x' = y, y' = -x - a*(x^2 - 1)*y
    J_vdp = [0, 1; -1, a_vdp];
    sys_vdp_origin = ss(J_vdp, zeros(2, 1), eye(2), zeros(2, 1));
    poles_vdp = pole(sys_vdp_origin);

    % Simulate van der Pol to convergence onto the unique Limit Cycle
    t_vdp = linspace(t_span_vdp(1), t_span_vdp(2), 3000);
    opts  = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    ode_vdp = @(t, U) [U(2); -U(1) - a_vdp * (U(1)^2 - 1.0) * U(2)];
    [~, U_vdp] = ode45(ode_vdp, t_vdp, [0.2; 0.0], opts);

    % Extract steady-state limit cycle metrics from the second half of the trajectory
    half_idx = round(length(t_vdp) / 2):length(t_vdp);
    t_ss = t_vdp(half_idx);
    x_ss = U_vdp(half_idx, 1);
    [pks, locs] = findpeaks(x_ss, t_ss);

    lc_amplitude = mean(pks);
    if length(locs) >= 2
        lc_period = mean(diff(locs));
    else
        lc_period = 2 * pi;
    end

    vdp_report.a_param = a_vdp;
    vdp_report.J_origin = J_vdp;
    vdp_report.origin_poles = poles_vdp;
    vdp_report.limit_cycle_amplitude = lc_amplitude;
    vdp_report.limit_cycle_period_s = lc_period;
    vdp_report.harmonic_period_2pi = 2 * pi;
    vdp_report.t = t_vdp;
    vdp_report.U_traj = U_vdp;

    %% 2. 3D Lorenz System Equilibria & Spectral Analysis (Section [4])
    a_lor = 10.0;
    b_lor = 8.0 / 3.0;

    % Origin Equilibrium E0 = (0, 0, 0)
    J_lor_origin = [-a_lor, a_lor,  0;
                     r_lorenz, -1,  0;
                     0,      0,    -b_lor];
    sys_lor_origin = ss(J_lor_origin, zeros(3, 1), eye(3), zeros(3, 1));
    poles_lor_origin = pole(sys_lor_origin);

    % Wing Equilibria C+/- (exist for r > 1)
    if r_lorenz > 1.0
        xy_eq = sqrt(b_lor * (r_lorenz - 1.0));
        z_eq  = r_lorenz - 1.0;
        J_lor_wing = [-a_lor,  a_lor,  0;
                       1.0,   -1.0,   -xy_eq;
                       xy_eq,  xy_eq, -b_lor];
        sys_lor_wing = ss(J_lor_wing, zeros(3, 1), eye(3), zeros(3, 1));
        poles_lor_wing = pole(sys_lor_wing);
    else
        xy_eq = 0; z_eq = 0;
        sys_lor_wing = sys_lor_origin;
        poles_lor_wing = poles_lor_origin;
    end

    t_lor = linspace(t_span_lor(1), t_span_lor(2), 4000);
    ode_lor = @(t, U) [-a_lor * U(1) + a_lor * U(2);
                       -U(1) * U(3) + r_lorenz * U(1) - U(2);
                        U(1) * U(2) - b_lor * U(3)];
    [~, U_lor] = ode45(ode_lor, t_lor, [1.0; 1.0; 20.0], opts);

    lorenz_report.parameters = struct('a', a_lor, 'b', b_lor, 'r', r_lorenz);
    lorenz_report.wing_equilibrium_pos = [xy_eq; xy_eq; z_eq];
    lorenz_report.origin_poles = poles_lor_origin;
    lorenz_report.wing_poles = poles_lor_wing;
    lorenz_report.is_wing_unstable = any(real(poles_lor_wing) > 0);
    lorenz_report.t = t_lor;
    lorenz_report.U_traj = U_lor;

    sys_models.sys_vdp_origin = sys_vdp_origin;
    sys_models.sys_lorenz_origin = sys_lor_origin;
    sys_models.sys_lorenz_wing = sys_lor_wing;

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Limit Cycles & Chaos Solver: MIT 18.03 Lec 39    \n');
    fprintf('====================================================\n');
    fprintf('--- [1] 2D van der Pol Oscillator (a = %.2f) ---\n', a_vdp);
    fprintf('Origin Linearized Poles  : %+.4f %+.4fi (Unstable Repeller)\n', ...
            real(poles_vdp(1)), abs(imag(poles_vdp(1))));
    fprintf('Steady Limit Cycle Amp   : X_max = %.4f (Theoretical ~ 2.0)\n', lc_amplitude);
    fprintf('Steady Limit Cycle Period: T_LC  = %.4f s (vs. 2*pi = %.4f s for a=0)\n', ...
            lc_period, 2*pi);
    fprintf('\n--- [2] 3D Lorenz System (r = %.2f) ---\n', r_lorenz);
    fprintf('Origin (0,0,0) Poles     : %+.3f, %+.3f, %+.3f\n', ...
            poles_lor_origin(1), poles_lor_origin(2), poles_lor_origin(3));
    fprintf('Wing Equilibrium C+      : (%.3f, %.3f, %.3f)\n', xy_eq, xy_eq, z_eq);
    fprintf('Wing C+/- Poles          : %+.3f, %+.3f %+.3fi (Unstable = %d)\n', ...
            real(poles_lor_wing(1)), real(poles_lor_wing(2)), abs(imag(poles_lor_wing(2))), ...
            lorenz_report.is_wing_unstable);
    fprintf('====================================================\n\n');

    %% 4. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [80, 80, 1180, 480]);

        % Left Subplot: van der Pol Limit Cycle & Time Waveform
        subplot(1, 2, 1);
        hold on; grid on; box on;
        plot(t_vdp, U_vdp(:, 1), 'b-', 'LineWidth', 1.8, 'DisplayName', 'Position x(t)');
        plot(t_vdp, U_vdp(:, 2), 'r--', 'LineWidth', 1.4, 'DisplayName', 'Velocity y(t) = x''(t)');
        yline(lc_amplitude, 'k:', 'LineWidth', 1.2, ...
              'DisplayName', sprintf('Limit Cycle Amp = %.2f', lc_amplitude));
        yline(-lc_amplitude, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('van der Pol Self-Sustained Oscillation (a=%.1f, T=%.2fs)', ...
              a_vdp, lc_period), 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'southeast', 'FontSize', 8.5);
        xlim([t_span_vdp(1), t_span_vdp(2)]);

        % Right Subplot: 3D Lorenz Strange Attractor
        subplot(1, 2, 2);
        plot3(U_lor(:, 1), U_lor(:, 2), U_lor(:, 3), 'b-', 'LineWidth', 1.0, ...
              'DisplayName', 'Strange Attractor Orbit');
        hold on; grid on; box on;
        if r_lorenz > 1.0
            plot3([xy_eq, -xy_eq], [xy_eq, -xy_eq], [z_eq, z_eq], 'ro', ...
                  'MarkerFaceColor', 'r', 'MarkerSize', 7, 'DisplayName', 'Wing Equilibria C_{\pm}');
        end
        xlabel('x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('y(t)', 'FontSize', 10, 'FontWeight', 'bold');
        zlabel('z(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('3D Lorenz Strange Attractor (r = %.1f)', r_lorenz), ...
              'FontSize', 11, 'FontWeight', 'bold');
        view([-38, 24]);
        legend('Location', 'northwest', 'FontSize', 8.5);
        hold off;
    end
end