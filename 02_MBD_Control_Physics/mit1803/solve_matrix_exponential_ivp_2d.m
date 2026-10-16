function [stm_meta, sys_models, sim_data] = solve_matrix_exponential_ivp_2d(A_mat, u0, Ts_sample, t_span)
% SOLVE_MATRIX_EXPONENTIAL_IVP_2D Computes the Matrix Exponential e^{At} = Phi(t)*Phi(0)^{-1}
% and solves the Initial Value Problem u' = A*u, u(0) = u0 across all 2x2 spectral cases.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (ss, initial, c2d, pole).
%
% % MBD Concept Mapping: Constructs the continuous State Transition Matrix (STM) e^{At}
% %                      via Fundamental Matrix normalization Phi(t)*Phi(0)^{-1}, and links it
% %                      to discrete-time ZOH state propagation u[k+1] = e^{A*Ts}*u[k] via 'c2d'.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 36 Section [2] Benchmark):
%   A_mat     - 2x2 system matrix A (default: [1.0, 1.0; 0.0, 2.0])
%   u0        - 2x1 initial state vector [x(0); y(0)] (default: [2.0; 1.0])
%   Ts_sample - Discrete ZOH sampling period for e^{A*Ts} evaluation (default: 0.25)
%   t_span    - [t_start, t_end] simulation interval (default: [0.0, 1.5])
%
% Outputs:
%   stm_meta   - Struct containing Phi(0), Phi(0)^{-1}, e^{A*Ts}, and analytical handlesPhi(t), e^{At}
%   sys_models - Struct containing continuous ss object and ZOH discretized ss object
%   sim_data   - Struct containing time vector, state trajectory u(t), and verification error

    arguments
        A_mat (2,2) double {mustBeNumeric} = [1.0, 1.0; 0.0, 2.0]
        u0 (2,1) double {mustBeNumeric} = [2.0; 1.0]
        Ts_sample (1,1) double {mustBePositive} = 0.25
        t_span (1,2) double {mustBeNumeric} = [0.0, 1.5]
    end

    %% 1. Control System Toolbox Continuous & Discrete State-Space Models
    sys_cont = ss(A_mat, zeros(2, 1), eye(2), zeros(2, 1));
    sys_disc = c2d(sys_cont, Ts_sample, 'zoh');

    sys_models.sys_continuous = sys_cont;
    sys_models.sys_discrete_zoh = sys_disc;

    %% 2. Construct Fundamental Matrix Phi(t) Across All Spectral Cases
    tr_A  = trace(A_mat);
    det_A = det(A_mat);
    disc  = tr_A^2 - 4 * det_A;
    tol   = 1e-11;

    if disc > tol
        % Distinct Real Eigenvalues: Phi(t) = [e^{l1*t}*v1, e^{l2*t}*v2]
        [V_raw, D_eig] = eig(A_mat);
        l1 = D_eig(1,1); l2 = D_eig(2,2);
        v1 = V_raw(:, 1); if abs(v1(1)) > 1e-9, v1 = v1 / v1(1); end
        v2 = V_raw(:, 2); if abs(v2(1)) > 1e-9, v2 = v2 / v2(1); end
        
        Phi_handle = @(t) [v1 * exp(l1 * t), v2 * exp(l2 * t)];
        case_str = "Distinct Real Eigenvalues";

    elseif disc < -tol
        % Complex Conjugate Eigenvalues: lambda = alpha +/- i*beta
        alpha = tr_A / 2.0;
        beta  = sqrt(-disc) / 2.0;
        M1 = A_mat - (alpha + 1i * beta) * eye(2);
        if norm(M1(1,:)) >= norm(M1(2,:))
            vc = [M1(1,2); -M1(1,1)];
        else
            vc = [M1(2,2); -M1(2,1)];
        end
        if abs(vc(1)) > 1e-9, vc = vc / vc(1); end
        vre = real(vc); vim = imag(vc);
        
        Phi_handle = @(t) exp(alpha * t) .* [vre*cos(beta*t) - vim*sin(beta*t), ...
                                             vre*sin(beta*t) + vim*cos(beta*t)];
        case_str = "Complex Conjugate Eigenvalues";

    else
        % Repeated Real Eigenvalue: lambda = tr(A) / 2
        l_rep = tr_A / 2.0;
        M_rep = A_mat - l_rep * eye(2);
        if norm(M_rep, 'fro') < tol
            % Complete case
            Phi_handle = @(t) exp(l_rep * t) * eye(2);
            case_str = "Repeated Eigenvalue (Complete Diagonal)";
        else
            % Defective case: u1 = e^{lt}*v, u2 = e^{lt}*(t*v + w)
            if norm(M_rep(1,:)) >= norm(M_rep(2,:))
                v_eig = [M_rep(1,2); -M_rep(1,1)];
            else
                v_eig = [M_rep(2,2); -M_rep(2,1)];
            end
            if abs(v_eig(1)) > 1e-9, v_eig = v_eig / v_eig(1); else, v_eig = v_eig / v_eig(2); end
            
            if abs(M_rep(1,2)) > 1e-9
                w_gen = [0; v_eig(1) / M_rep(1,2)];
            else
                w_gen = pinv(M_rep) * v_eig;
            end
            Phi_handle = @(t) exp(l_rep * t) .* [v_eig, v_eig * t + w_gen];
            case_str = "Repeated Eigenvalue (Defective Jordan Block)";
        end
    end

    %% 3. Normalize to Obtain Matrix Exponential: e^{At} = Phi(t) * Phi(0)^{-1}
    Phi_0 = Phi_handle(0);
    det_Phi_0 = det(Phi_0);
    % Explicit 2x2 inverse formula from Lecture Section [2]
    Phi_0_inv = (1.0 / det_Phi_0) * [ Phi_0(2,2), -Phi_0(1,2); 
                                     -Phi_0(2,1),  Phi_0(1,1)];

    eAt_handle = @(t) Phi_handle(t) * Phi_0_inv;
    u_sol_handle = @(t) eAt_handle(t) * u0;

    stm_meta.spectral_case = case_str;
    stm_meta.Phi_0 = Phi_0;
    stm_meta.Phi_0_inv = Phi_0_inv;
    stm_meta.Phi_handle = Phi_handle;
    stm_meta.eAt_handle = eAt_handle;
    stm_meta.u_sol_handle = u_sol_handle;
    stm_meta.eA_at_t1 = eAt_handle(1.0);
    stm_meta.Ad_zoh = sys_disc.A;

    %% 4. Numerical Simulation & Column-by-Column Basis Verification
    t_eval = linspace(t_span(1), t_span(2), 500);
    [u_sim_cont, ~] = initial(sys_cont, u0, t_eval);

    % Evaluate analytical u(t) = e^{At} * u0 across time grid
    u_analytic = zeros(length(t_eval), 2);
    eAt_entries = zeros(length(t_eval), 4); % [e11, e12, e21, e22]
    for idx = 1:length(t_eval)
        M_exp = eAt_handle(t_eval(idx));
        u_analytic(idx, :) = (M_exp * u0).';
        eAt_entries(idx, :) = [M_exp(1,1), M_exp(1,2), M_exp(2,1), M_exp(2,2)];
    end

    % Discrete-time trajectory via Power Law: u[k] = (e^{A*Ts})^k * u0
    t_disc = t_span(1):Ts_sample:t_span(2);
    [u_sim_disc, ~] = initial(sys_disc, u0, t_disc);

    max_err = max(max(abs(u_sim_cont - u_analytic)));

    sim_data.t = t_eval;
    sim_data.u_sim = u_sim_cont;
    sim_data.u_analytic = u_analytic;
    sim_data.eAt_entries = eAt_entries;
    sim_data.t_disc = t_disc;
    sim_data.u_disc = u_sim_disc;
    sim_data.verification_max_error = max_err;

    %% 5. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Matrix Exponential Solver: MIT 18.03 Lecture 36  \n');
    fprintf('====================================================\n');
    fprintf('Spectral Case           : %s\n', case_str);
    fprintf('Fundamental Matrix Phi(0):\n');
    disp(Phi_0);
    fprintf('Inverse Matrix Phi(0)^{-1}:\n');
    disp(Phi_0_inv);
    fprintf('Matrix Exponential e^{A} (at t = 1s):\n');
    disp(stm_meta.eA_at_t1);
    fprintf('Verification Max Error  : %e (|ss initial - Phi(t)Phi(0)^{-1}u0|)\n', max_err);
    fprintf('====================================================\n\n');

    %% 6. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [100, 100, 1120, 460]);

        % Left Subplot: Evolution of the 4 Matrix Exponential Entries [e^{At}]_{ij}
        subplot(1, 2, 1);
        hold on; grid on; box on;
        plot(t_eval, eAt_entries(:, 1), 'b-',  'LineWidth', 1.8, 'DisplayName', '[e^{At}]_{11} (Col 1, Row 1)');
        plot(t_eval, eAt_entries(:, 2), 'r-',  'LineWidth', 1.8, 'DisplayName', '[e^{At}]_{12} (Col 2, Row 1)');
        plot(t_eval, eAt_entries(:, 3), 'g--', 'LineWidth', 1.8, 'DisplayName', '[e^{At}]_{21} (Col 1, Row 2)');
        plot(t_eval, eAt_entries(:, 4), 'm-.', 'LineWidth', 1.8, 'DisplayName', '[e^{At}]_{22} (Col 2, Row 2)');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Matrix Entry Value', 'FontSize', 10, 'FontWeight', 'bold');
        title('State Transition Matrix Entries: e^{At} = \Phi(t)\Phi(0)^{-1}', ...
              'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);

        % Right Subplot: Continuous vs. Discrete Power-Law State Propagation u(t) = e^{At}u(0)
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, u_sim_cont(:, 1), 'b-', 'LineWidth', 2.0, 'DisplayName', 'Continuous x(t) = [e^{At}u_0]_1');
        plot(t_eval, u_sim_cont(:, 2), 'r-', 'LineWidth', 2.0, 'DisplayName', 'Continuous y(t) = [e^{At}u_0]_2');
        stairs(t_disc, u_sim_disc(:, 1), 'k--', 'LineWidth', 1.2, ...
               'DisplayName', sprintf('Discrete ZOH (e^{AT_s})^k, T_s=%.2fs', Ts_sample));
        plot(t_disc, u_sim_disc(:, 1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
        plot(t_disc, u_sim_disc(:, 2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'HandleVisibility', 'off');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Amplitude u(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('IVP Solution from u(0) = [%.1f; %.1f]', u0(1), u0(2)), ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end