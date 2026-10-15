function [sys_meta, sys_ss, sim_data] = solve_linear_first_order_system_2d(A_or_coeffs, u0, t_span)
% SOLVE_LINEAR_FIRST_ORDER_SYSTEM_2D Solves 2D linear systems u' = A*u and Companion Matrices.
% Supports both direct 2x2 matrix A and 2nd-order scalar ODE coefficients [1, b, k].
%
% % Required Toolboxes: 
% %   - Control System Toolbox (ss, initial, pole).
%
% % MBD Concept Mapping: Bridges 2x2 State-Space systems (u' = A*u) and 2nd-order scalar ODEs
% %                      via Elimination (ss -> tf) and Anti-Elimination Companion Matrices (tf -> ss).
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 32 Romeo & Juliet Benchmark):
%   A_or_coeffs - Either a 2x2 matrix A [a, b; c, d] OR a 1x3 vector [m, b, k] for m*x'' + b*x' + k*x = 0
%                 (default: [0, 1; -17/16, 1/2])
%   u0          - Initial state column vector [x(0); y(0)] (default: [1.0; 0.0])
%   t_span      - [t_start, t_end] simulation interval (default: [0, 4*pi])
%
% Outputs:
%   sys_meta - Struct with matrix A, trace, det, characteristic poly, poles, eigenvectors, and topology
%   sys_ss   - State-Space object from Control System Toolbox
%   sim_data - Struct with time vector, state trajectories [x(t), y(t)], and verification error

    arguments
        A_or_coeffs double {mustBeNumeric} = [0, 1; -17/16, 1/2]
        u0 (2,1) double {mustBeNumeric} = [1.0; 0.0]
        t_span (1,2) double {mustBeNumeric} = [0, 4 * pi]
    end

    %% 1. Parse Input: Direct 2x2 Matrix vs. Companion Matrix Construction
    if isequal(size(A_or_coeffs), [2, 2])
        A_mat = A_or_coeffs;
        sys_meta.input_mode = "Direct 2x2 System Matrix A";
    elseif numel(A_or_coeffs) == 3
        % Anti-elimination: m*x'' + b*x' + k*x = 0 -> Companion matrix [0, 1; -k/m, -b/m]
        m_val = A_or_coeffs(1); b_val = A_or_coeffs(2); k_val = A_or_coeffs(3);
        A_mat = [0, 1; -k_val / m_val, -b_val / m_val];
        sys_meta.input_mode = "Companion Matrix from 2nd-Order ODE";
    else
        error('Input must be a 2x2 matrix A or a 3-element ODE coefficient vector [m, b, k].');
    end

    %% 2. State-Space Realization & Eigen-Decomposition
    sys_ss = ss(A_mat, zeros(2, 1), eye(2), zeros(2, 1));
    
    tr_A  = trace(A_mat);
    det_A = det(A_mat);
    char_poly = [1.0, -tr_A, det_A]; % Eliminated 2nd-order ODE: x'' - tr(A)*x' + det(A)*x = 0
    
    [V_eig, D_eig] = eig(A_mat);
    eigenvals = diag(D_eig);

    % Classify Phase Portrait Topology
    disc = tr_A^2 - 4 * det_A;
    if disc < -1e-12
        if abs(tr_A) < 1e-9
            topology_str = "Center (Closed Elliptical Orbits)";
        elseif tr_A > 0
            topology_str = "Unstable Spiral / Focus (Outward Expanding Oscillations)";
        else
            topology_str = "Stable Spiral / Focus (Inward Decaying Oscillations)";
        end
    elseif det_A < 0
        topology_str = "Saddle Point (One Stable Mode, One Unstable Mode)";
    else
        if tr_A > 0
            topology_str = "Unstable Node (Exponential Growth along Eigenvectors)";
        else
            topology_str = "Stable Node (Exponential Decay along Eigenvectors)";
        end
    end

    sys_meta.A = A_mat;
    sys_meta.trace = tr_A;
    sys_meta.determinant = det_A;
    sys_meta.characteristic_poly = char_poly;
    sys_meta.eigenvalues = eigenvals;
    sys_meta.eigenvectors = V_eig;
    sys_meta.topology = topology_str;

    %% 3. Time-Domain Simulation via Control System Toolbox 'initial' & Matrix Exponential Check
    t_eval = linspace(t_span(1), t_span(2), 1000);
    [u_sim, ~] = initial(sys_ss, u0, t_eval);

    % Exact analytical matrix exponential evaluation: u(t) = expm(A*t) * u0
    u_exact = zeros(length(t_eval), 2);
    for idx = 1:length(t_eval)
        u_exact(idx, :) = (expm(A_mat * t_eval(idx)) * u0).';
    end

    max_err = max(max(abs(u_sim - u_exact)));

    sim_data.t = t_eval;
    sim_data.u_sim = u_sim;
    sim_data.u_exact = u_exact;
    sim_data.verification_max_error = max_err;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   2D Linear System Solver: MIT 18.03 Lecture 32    \n');
    fprintf('====================================================\n');
    fprintf('System Mode             : %s\n', sys_meta.input_mode);
    fprintf('Matrix A                : [ %+.4f, %+.4f ; %+.4f, %+.4f ]\n', ...
            A_mat(1,1), A_mat(1,2), A_mat(2,1), A_mat(2,2));
    fprintf('Eliminated 2nd-Order ODE: x'''' %+.4f*x'' %+.4f*x = 0\n', -tr_A, det_A);
    fprintf('Eigenvalues (Poles)     : r1 = %+.4f %+.4fi,  r2 = %+.4f %+.4fi\n', ...
            real(eigenvals(1)), imag(eigenvals(1)), real(eigenvals(2)), imag(eigenvals(2)));
    fprintf('Phase Plane Topology    : %s\n', topology_str);
    fprintf('Verification Max Error  : %e (|ss initial - expm(At)u0|)\n', max_err);
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [100, 100, 1120, 460]);

        % Left Subplot: Phase Plane Trajectory & Vector Field
        subplot(1, 2, 1);
        hold on; grid on; box on;
        xline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        max_bound = max(max(abs(u_sim))) * 1.1;
        if max_bound == 0, max_bound = 1; end
        [X_grid, Y_grid] = meshgrid(linspace(-max_bound, max_bound, 15));
        U_vel = A_mat(1,1) * X_grid + A_mat(1,2) * Y_grid;
        V_vel = A_mat(2,1) * X_grid + A_mat(2,2) * Y_grid;
        quiver(X_grid, Y_grid, U_vel, V_vel, 'Color', [0.75 0.75 0.75], ...
               'LineWidth', 0.8, 'DisplayName', 'Velocity Field A*u');

        % If eigenvalues are real, plot invariant eigenvector lines
        if disc >= 0
            s_ray = linspace(-max_bound, max_bound, 100);
            plot(s_ray * real(V_eig(1,1)), s_ray * real(V_eig(2,1)), 'r--', ...
                 'LineWidth', 1.4, 'DisplayName', sprintf('Mode 1 Ray (r=%.2f)', real(eigenvals(1))));
            plot(s_ray * real(V_eig(1,2)), s_ray * real(V_eig(2,2)), 'g--', ...
                 'LineWidth', 1.4, 'DisplayName', sprintf('Mode 2 Ray (r=%.2f)', real(eigenvals(2))));
        end

        plot(u_sim(:, 1), u_sim(:, 2), 'b-', 'LineWidth', 2.0, 'DisplayName', 'Trajectory u(t) = [x(t); y(t)]');
        plot(u0(1), u0(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'DisplayName', 'Initial State u(0)');
        plot(u_sim(end, 1), u_sim(end, 2), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'DisplayName', 'Final State');

        xlabel('State x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State y(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Phase Plane: %s', topology_str), 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8);
        xlim([-max_bound, max_bound]); ylim([-max_bound, max_bound]);

        % Right Subplot: Time-Domain State Evolution
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, u_sim(:, 1), 'b-', 'LineWidth', 2.0, 'DisplayName', 'State x(t) (ss initial)');
        plot(t_eval, u_sim(:, 2), 'r-', 'LineWidth', 1.8, 'DisplayName', 'State y(t) (ss initial)');
        plot(t_eval(1:40:end), u_exact(1:40:end, 1), 'ko', 'MarkerSize', 5, 'DisplayName', 'Matrix Exp expm(At)u_0');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title('Time-Domain Evolution of Coupled States', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end