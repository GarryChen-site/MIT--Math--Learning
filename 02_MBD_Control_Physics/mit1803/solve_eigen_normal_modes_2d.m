function [modal_decomp, sys_ss, sim_data] = solve_eigen_normal_modes_2d(A_mat, u0, t_span)
% SOLVE_EIGEN_NORMAL_MODES_2D Solves u' = A*u via Eigenvalues, Eigenvectors & Normal Modes.
% Decomposes u(t) = c1*e^{lambda1*t}*v1 + c2*e^{lambda2*t}*v2 and validates via Control Toolbox.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (ss, initial, pole).
%
% % MBD Concept Mapping: Extracts system poles (eigenvalues) and mode shapes (eigenvectors).
% %                      Performs modal coordinate transformation c = V \ u0 and identifies
% %                      asymptotic dominant manifolds as t -> +infty and t -> -infty.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 33 Benchmark):
%   A_mat  - 2x2 system matrix A (default: [-2, 1; -4, 3])
%   u0     - 2x1 initial state column vector [x(0); y(0)] (default: [1.0; 0.0])
%   t_span - [t_start, t_end] forward simulation interval (default: [0, 1.5])
%
% Outputs:
%   modal_decomp - Struct containing eigenvalues, normalized eigenvectors, c1, c2, and handles
%   sys_ss       - Autonomous State-Space object from Control System Toolbox
%   sim_data     - Struct containing time vector, modal trajectories, and verification error

    arguments
        A_mat (2,2) double {mustBeNumeric} = [-2.0, 1.0; -4.0, 3.0]
        u0 (2,1) double {mustBeNumeric} = [1.0; 0.0]
        t_span (1,2) double {mustBeNumeric} = [0.0, 1.5]
    end

    %% 1. State-Space Model & Characteristic Polynomial
    sys_ss = ss(A_mat, zeros(2, 1), eye(2), zeros(2, 1));
    
    tr_A  = trace(A_mat);
    det_A = det(A_mat);
    char_poly = [1.0, -tr_A, det_A]; % p_A(lambda) = lambda^2 - tr(A)*lambda + det(A)

    %% 2. Eigenvalue & Clean Eigenvector Extraction
    [V_raw, D_eig] = eig(A_mat);
    lambda_vec = diag(D_eig);

    % Normalize eigenvectors so that the first non-zero component is 1 (matches hand calculations)
    V_clean = V_raw;
    for k = 1:2
        if abs(V_raw(1, k)) > 1e-9
            V_clean(:, k) = V_raw(:, k) / V_raw(1, k);
        elseif abs(V_raw(2, k)) > 1e-9
            V_clean(:, k) = V_raw(:, k) / V_raw(2, k);
        end
    end

    % Solve for modal superposition weights: V_clean * [c1; c2] = u0
    if rank(V_clean) == 2
        c_weights = V_clean \ u0;
    else
        warning('Defective matrix (repeated eigenvalue with single eigenvector). Using raw matrix exponential.');
        c_weights = [NaN; NaN];
    end

    modal_decomp.A = A_mat;
    modal_decomp.characteristic_poly = char_poly;
    modal_decomp.eigenvalues = lambda_vec;
    modal_decomp.eigenvectors = V_clean;
    modal_decomp.c_weights = c_weights;

    % Analytical function handles for the two normal modes and total solution
    v1 = V_clean(:, 1); l1 = lambda_vec(1); c1 = c_weights(1);
    v2 = V_clean(:, 2); l2 = lambda_vec(2); c2 = c_weights(2);

    modal_decomp.mode1_handle = @(t) c1 * exp(l1 * t) .* v1;
    modal_decomp.mode2_handle = @(t) c2 * exp(l2 * t) .* v2;
    modal_decomp.u_total_handle = @(t) real(c1 * exp(l1 * t) .* v1 + c2 * exp(l2 * t) .* v2);

    %% 3. Forward Time Simulation via Control System Toolbox 'initial'
    t_eval = linspace(t_span(1), t_span(2), 600);
    [u_sim, ~] = initial(sys_ss, u0, t_eval);

    if ~any(isnan(c_weights))
        u_analytic = modal_decomp.u_total_handle(t_eval).';
    else
        u_analytic = zeros(length(t_eval), 2);
        for idx = 1:length(t_eval)
            u_analytic(idx, :) = (expm(A_mat * t_eval(idx)) * u0).';
        end
    end

    max_err = max(max(abs(u_sim - u_analytic)));

    sim_data.t = t_eval;
    sim_data.u_sim = u_sim;
    sim_data.u_analytic = u_analytic;
    sim_data.verification_max_error = max_err;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Eigenvalue & Normal Mode Solver: MIT 18.03 Lec 33\n');
    fprintf('====================================================\n');
    fprintf('Characteristic Poly p(lambda): lambda^2 %+.4f*lambda %+.4f = 0\n', -tr_A, det_A);
    for k = 1:2
        fprintf('Mode %d: lambda_%d = %+.4f %+.4fi | v_%d = [%+.4f; %+.4f] | c_%d = %+.4f\n', ...
                k, k, real(lambda_vec(k)), imag(lambda_vec(k)), ...
                k, real(V_clean(1,k)), real(V_clean(2,k)), k, real(c_weights(k)));
    end
    fprintf('Verification Max Error       : %e (|ss initial - modal sum|)\n', max_err);
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [100, 100, 1120, 460]);

        % Left Subplot: Phase Plane with Normal Mode Rays
        subplot(1, 2, 1);
        hold on; grid on; box on;
        xline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        max_val = max(max(abs(u_sim))) * 1.1;
        if max_val == 0, max_val = 1; end

        if isreal(V_clean)
            s_line = linspace(-max_val, max_val, 100);
            plot(s_line * v1(1), s_line * v1(2), 'b--', 'LineWidth', 1.6, ...
                 'DisplayName', sprintf('Normal Mode 1 Ray ($\\lambda_1=%.2f$)', l1));
            plot(s_line * v2(1), s_line * v2(2), 'r--', 'LineWidth', 1.6, ...
                 'DisplayName', sprintf('Normal Mode 2 Ray ($\\lambda_2=%.2f$)', l2));
        end

        plot(u_sim(:, 1), u_sim(:, 2), 'k-', 'LineWidth', 2.2, 'DisplayName', 'Trajectory u(t) ($t \ge 0$)');
        plot(u0(1), u0(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Initial State u(0)');
        plot(u_sim(end, 1), u_sim(end, 2), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'DisplayName', 'Final State');

        xlabel('State x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State y(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title('Phase Plane: Normal Mode Rays & Trajectory', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5, 'Interpreter','latex');
        xlim([-max_val, max_val]); ylim([-max_val, max_val]);

        % Right Subplot: Time-Domain State Trajectories
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, u_sim(:, 1), 'b-', 'LineWidth', 2.0, 'DisplayName', 'State x(t) (ss initial)');
        plot(t_eval, u_sim(:, 2), 'r-', 'LineWidth', 2.0, 'DisplayName', 'State y(t) (ss initial)');
        plot(t_eval(1:30:end), u_analytic(1:30:end, 1), 'ko', 'MarkerSize', 5, 'DisplayName', 'Modal Sum x(t)');
        plot(t_eval(1:30:end), u_analytic(1:30:end, 2), 'go', 'MarkerSize', 5, 'DisplayName', 'Modal Sum y(t)');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title('Time Domain: Modal Superposition Verification', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end