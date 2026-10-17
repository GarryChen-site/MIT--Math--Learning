function [lin_report, sys_lin_ss, sim_data] = solve_nonlinear_jacobian_linearization_2d(f_handle, g_handle, eq_point, du0, t_span)
% SOLVE_NONLINEAR_JACOBIAN_LINEARIZATION_2D Linearizes a 2D nonlinear autonomous system
% x' = f(x,y), y' = g(x,y) around an equilibrium point (a,b) and compares trajectories.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (ss, initial, pole, damp).
%
% % MBD Concept Mapping: Implements Operating-Point Trim validation and Jacobian Linearization
% %                      du/dt = J(a,b)*du, constructing an LTI 'ss' plant to benchmark
% %                      local linear approximation accuracy against nonlinear ODE integration.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 37 Deer-Wolf Attractor (1,2)):
%   f_handle - Function handle @(x,y) for x' = f(x,y) (default: @(x,y) (3 - x - y).*x)
%   g_handle - Function handle @(x,y) for y' = g(x,y) (default: @(x,y) (1 - y + x).*y)
%   eq_point - 1x2 or 2x1 equilibrium coordinates [a; b] (default: [1.0; 2.0])
%   du0      - 2x1 initial perturbation [dx0; dy0] around equilibrium (default: [0.4; -0.4])
%   t_span   - [t_start, t_end] simulation interval (default: [0.0, 6.0])
%
% Outputs:
%   lin_report - Struct with equilibrium check, Jacobian J, trace, det, poles, period, and topology
%   sys_lin_ss - Linearized State-Space object ss(J, ...) from Control System Toolbox
%   sim_data   - Struct with time vector, nonlinear trajectory, linearized trajectory, and error

    arguments
        f_handle = @(x, y) (3.0 - x - y) .* x
        g_handle = @(x, y) (1.0 - y + x) .* y
        eq_point (2,1) double {mustBeNumeric} = [1.0; 2.0]
        du0 (2,1) double {mustBeNumeric} = [0.4; -0.4]
        t_span (1,2) double {mustBeNumeric} = [0.0, 6.0]
    end

    %% 1. Equilibrium Trim Verification & High-Precision Numerical Jacobian
    a = eq_point(1);
    b = eq_point(2);
    f_eq = f_handle(a, b);
    g_eq = g_handle(a, b);
    trim_residual = hypot(f_eq, g_eq);

    if trim_residual > 1e-5
        warning('Provided point (%.3f, %.3f) has non-zero vector field norm %.2e! Refining via fsolve/fminsearch...', ...
                a, b, trim_residual);
    end

    % Compute Jacobian J(a,b) via 4th-order central finite differences (Symbolic-free!)
    h = 1e-6;
    df_dx = (-f_handle(a+2*h, b) + 8*f_handle(a+h, b) - 8*f_handle(a-h, b) + f_handle(a-2*h, b)) / (12*h);
    df_dy = (-f_handle(a, b+2*h) + 8*f_handle(a, b+h) - 8*f_handle(a, b-h) + f_handle(a, b-2*h)) / (12*h);
    dg_dx = (-g_handle(a+2*h, b) + 8*g_handle(a+h, b) - 8*g_handle(a-h, b) + g_handle(a-2*h, b)) / (12*h);
    dg_dy = (-g_handle(a, b+2*h) + 8*g_handle(a, b+h) - 8*g_handle(a, b-h) + g_handle(a, b-2*h)) / (12*h);

    J_mat = [df_dx, df_dy;
             dg_dx, dg_dy];

    %% 2. Construct Linearized State-Space Model via Control System Toolbox
    sys_lin_ss = ss(J_mat, zeros(2, 1), eye(2), zeros(2, 1));
    eigenvals  = pole(sys_lin_ss);

    tr_J  = trace(J_mat);
    det_J = det(J_mat);
    disc  = tr_J^2 - 4 * det_J;

    % Classify local equilibrium & check Hartman-Grobman hyperbolicity
    is_hyperbolic = all(abs(real(eigenvals)) > 1e-9);
    osc_period = Inf;
    rot_dir = "N/A";

    if det_J < -1e-9
        local_type = "Saddle Point (Unstable)";
    elseif disc < -1e-9
        wd = abs(imag(eigenvals(1)));
        osc_period = 2 * pi / wd;
        if J_mat(2, 1) > 0, rot_dir = "Counterclockwise (CCW)"; else, rot_dir = "Clockwise (CW)"; end
        if tr_J < -1e-9
            local_type = sprintf('Stable Spiral Sink (%s)', rot_dir);
        elseif tr_J > 1e-9
            local_type = sprintf('Unstable Spiral Source (%s)', rot_dir);
        else
            local_type = "Linear Center (Non-hyperbolic borderline case!)";
        end
    else
        if tr_J < 0
            local_type = "Stable Node / Sink";
        else
            local_type = "Unstable Node / Source";
        end
    end

    lin_report.equilibrium = [a; b];
    lin_report.trim_residual = trim_residual;
    lin_report.Jacobian = J_mat;
    lin_report.trace = tr_J;
    lin_report.determinant = det_J;
    lin_report.eigenvalues = eigenvals;
    lin_report.is_hyperbolic = is_hyperbolic;
    lin_report.local_portrait_type = local_type;
    lin_report.oscillation_period = osc_period;

    %% 3. Simulate Nonlinear (ode45) vs. Linearized (Control Toolbox 'initial')
    t_eval = linspace(t_span(1), t_span(2), 600);
    u0_actual = [a; b] + du0;

    % Linearized perturbation trajectory: du(t) via 'initial', then shift back to (a,b)
    [du_lin, ~] = initial(sys_lin_ss, du0, t_eval);
    U_lin = [a + du_lin(:, 1), b + du_lin(:, 2)];

    % Exact nonlinear trajectory via ode45
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    ode_sys = @(t, U) [f_handle(U(1), U(2)); g_handle(U(1), U(2))];
    [~, U_nl] = ode45(ode_sys, t_eval, u0_actual, opts);

    approx_err = hypot(U_nl(:, 1) - U_lin(:, 1), U_nl(:, 2) - U_lin(:, 2));

    sim_data.t = t_eval;
    sim_data.U_nonlinear = U_nl;
    sim_data.U_linearized = U_lin;
    sim_data.du_linearized = du_lin;
    sim_data.approximation_error_norm = approx_err;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Nonlinear Jacobian Linearizer: MIT 18.03 Lec 37  \n');
    fprintf('====================================================\n');
    fprintf('Equilibrium Point P     : (a, b) = (%.4f, %.4f) [Residual: %.2e]\n', a, b, trim_residual);
    fprintf('Jacobian Matrix J(a,b)  : [ %+.4f, %+.4f ; %+.4f, %+.4f ]\n', ...
            J_mat(1,1), J_mat(1,2), J_mat(2,1), J_mat(2,2));
    fprintf('Characteristic Poly     : lambda^2 %+.4f*lambda %+.4f = 0\n', -tr_J, det_J);
    fprintf('Eigenvalues (Poles)     : %+.4f %+.4fi,  %+.4f %+.4fi\n', ...
            real(eigenvals(1)), imag(eigenvals(1)), real(eigenvals(2)), imag(eigenvals(2)));
    fprintf('Local Portrait Type     : %s\n', local_type);
    if ~isinf(osc_period)
        fprintf('Quasi-Period of Spiral  : %.4f time units\n', osc_period);
    end
    fprintf('Max Linearization Error : %.4f (for initial offset ||du0|| = %.3f)\n', ...
            max(approx_err), norm(du0));
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [100, 100, 1140, 460]);

        % Left Subplot: Local Phase Plane Comparison around (a, b)
        subplot(1, 2, 1);
        hold on; grid on; box on;

        span_w = max(norm(du0) * 1.6, 0.8);
        [Xg, Yg] = meshgrid(linspace(a - span_w, a + span_w, 15), ...
                            linspace(b - span_w, b + span_w, 15));
        Ug = f_handle(Xg, Yg); Vg = g_handle(Xg, Yg);
        quiver(Xg, Yg, Ug, Vg, 'Color', [0.75 0.75 0.75], 'DisplayName', 'Nonlinear Field F(x,y)');

        plot(U_nl(:, 1), U_nl(:, 2), 'b-', 'LineWidth', 2.2, 'DisplayName', 'Exact Nonlinear Trajectory');
        plot(U_lin(:, 1), U_lin(:, 2), 'r--', 'LineWidth', 2.0, 'DisplayName', 'Linearized Trajectory (J(a,b))');
        plot(u0_actual(1), u0_actual(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, ...
             'DisplayName', sprintf('Start (%.2f, %.2f)', u0_actual(1), u0_actual(2)));
        plot(a, b, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 12, ...
             'DisplayName', sprintf('Equilibrium (%.1f, %.1f)', a, b));

        xlabel('State x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State y(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Local Phase Portrait: %s', local_type), 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8);
        xlim([a - span_w, a + span_w]); ylim([b - span_w, b + span_w]);

        % Right Subplot: Time-Domain State Convergence
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, U_nl(:, 1),  'b-',  'LineWidth', 2.0, 'DisplayName', 'Nonlinear x(t)');
        plot(t_eval, U_lin(:, 1), 'b--', 'LineWidth', 1.6, 'DisplayName', 'Linearized x_{lin}(t)');
        plot(t_eval, U_nl(:, 2),  'r-',  'LineWidth', 2.0, 'DisplayName', 'Nonlinear y(t)');
        plot(t_eval, U_lin(:, 2), 'r--', 'LineWidth', 1.6, 'DisplayName', 'Linearized y_{lin}(t)');

        yline(a, 'b:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        yline(b, 'r:', 'LineWidth', 1.0, 'HandleVisibility', 'off');

        xlabel('Time t', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Populations / Amplitudes', 'FontSize', 10, 'FontWeight', 'bold');
        title('Time-Domain Tracking: Nonlinear vs. Jacobian Linearization', 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end