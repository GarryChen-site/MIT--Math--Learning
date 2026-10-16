function [spectral_info, sys_ss, sim_data] = solve_complex_or_repeated_eigen_2d(A_mat, u0, t_span)
% SOLVE_COMPLEX_OR_REPEATED_EIGEN_2D Universal 2x2 solver for u' = A*u across all spectral cases:
%   1) Complex Conjugate Eigenvalues -> Real & Imaginary spiral basis u1 = Re(u), u2 = Im(u).
%   2) Defective Repeated Eigenvalue -> Generalized eigenvector (A - lambda*I)*w = v, u2 = e^{lt}*(t*v + w).
%   3) Complete Repeated Eigenvalue  -> Star node basis u1 = e^{lt}*[1;0], u2 = e^{lt}*[0;1].
%   4) Distinct Real Eigenvalues     -> Normal mode rays u1 = e^{l1*t}*v1, u2 = e^{l2*t}*v2.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (ss, initial, pole).
%
% % MBD Concept Mapping: Handles Jordan blocks (defective matrices) and Real Modal blocks (complex poles)
% %                      to guarantee a complete 2D basis for any initial state u(0).
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 34 Sec [3] Defective Benchmark):
%   A_mat  - 2x2 system matrix A (default: [-2, 1; -1, 0])
%   u0     - 2x1 initial state column vector [x(0); y(0)] (default: [0.0; 1.0])
%   t_span - [t_start, t_end] simulation interval (default: [0.0, 5.0])
%
% Outputs:
%   spectral_info - Struct with spectral case, eigenvalues, basis vectors (v, w), weights, and handles
%   sys_ss        - State-Space object from Control System Toolbox
%   sim_data      - Struct with time vector, simulated trajectory, analytical trajectory, and max error

    arguments
        A_mat (2,2) double {mustBeNumeric} = [-2.0, 1.0; -1.0, 0.0]
        u0 (2,1) double {mustBeNumeric} = [0.0; 1.0]
        t_span (1,2) double {mustBeNumeric} = [0.0, 5.0]
    end

    %% 1. Invariants & State-Space Construction via Control System Toolbox
    sys_ss = ss(A_mat, zeros(2, 1), eye(2), zeros(2, 1));
    
    tr_A  = trace(A_mat);
    det_A = det(A_mat);
    disc  = tr_A^2 - 4 * det_A;

    spectral_info.A = A_mat;
    spectral_info.trace = tr_A;
    spectral_info.determinant = det_A;
    spectral_info.characteristic_poly = [1.0, -tr_A, det_A];

    %% 2. Spectral Branching & Analytical Basis Construction
    if disc < -1e-11
        % =================================================================
        % Case 1: Complex Conjugate Eigenvalues (lambda = alpha +/- i*beta)
        % =================================================================
        alpha = tr_A / 2.0;
        beta  = sqrt(-disc) / 2.0;
        lambda1 = alpha + 1i * beta;
        
        % Construct complex eigenvector v1 for lambda1 from top or bottom row of (A - lambda1*I)
        M1 = A_mat - lambda1 * eye(2);
        if norm(M1(1, :)) >= norm(M1(2, :))
            v_complex = [M1(1, 2); -M1(1, 1)];
        else
            v_complex = [M1(2, 2); -M1(2, 1)];
        end
        % Normalize first component to 1 if non-zero (matches lecture v1 = [1; i])
        if abs(v_complex(1)) > 1e-9
            v_complex = v_complex / v_complex(1);
        end
        
        v_re = real(v_complex);
        v_im = imag(v_complex);
        
        % Real and Imaginary part basis solutions: u(t) = e^{(alpha + i*beta)t} * (v_re + i*v_im)
        u1_handle = @(t) exp(alpha * t) .* (cos(beta * t) .* v_re - sin(beta * t) .* v_im);
        u2_handle = @(t) exp(alpha * t) .* (sin(beta * t) .* v_re + cos(beta * t) .* v_im);
        
        spectral_info.case_type = "Complex Conjugate Eigenvalues (Spiral)";
        spectral_info.eigenvalues = [lambda1; conj(lambda1)];
        spectral_info.v_primary = v_complex;
        spectral_info.w_generalized = [NaN; NaN];

    elseif abs(disc) <= 1e-11
        % =================================================================
        % Repeated Eigenvalue: lambda1 = lambda2 = tr(A) / 2
        % =================================================================
        lambda_rep = tr_A / 2.0;
        M_rep = A_mat - lambda_rep * eye(2);
        
        if norm(M_rep, 'fro') < 1e-11
            % Case 2A: Complete Repeated Matrix A = [lambda, 0; 0, lambda]
            v1 = [1; 0];
            v2 = [0; 1];
            u1_handle = @(t) exp(lambda_rep * t) .* v1;
            u2_handle = @(t) exp(lambda_rep * t) .* v2;
            
            spectral_info.case_type = "Repeated Eigenvalue - Complete (Star Node)";
            spectral_info.eigenvalues = [lambda_rep; lambda_rep];
            spectral_info.v_primary = v1;
            spectral_info.w_generalized = v2;
        else
            % Case 2B: Defective Repeated Matrix (Only 1 Eigenline)
            if norm(M_rep(1, :)) >= norm(M_rep(2, :))
                v_eig = [M_rep(1, 2); -M_rep(1, 1)];
            else
                v_eig = [M_rep(2, 2); -M_rep(2, 1)];
            end
            if abs(v_eig(1)) > 1e-9
                v_eig = v_eig / v_eig(1); % Normalize to match [1; 1]
            else
                v_eig = v_eig / v_eig(2);
            end
            
            % Find generalized eigenvector w such that (A - lambda*I)*w = v_eig
            % Match lecture convention: try setting w(1) = 0 first if M_rep(1,2) ~= 0
            if abs(M_rep(1, 2)) > 1e-9
                w_gen = [0; v_eig(1) / M_rep(1, 2)];
            else
                w_gen = pinv(M_rep) * v_eig;
            end
            
            u1_handle = @(t) exp(lambda_rep * t) .* v_eig;
            u2_handle = @(t) exp(lambda_rep * t) .* (v_eig * t + w_gen);
            
            spectral_info.case_type = "Repeated Eigenvalue - Defective (Improper Node)";
            spectral_info.eigenvalues = [lambda_rep; lambda_rep];
            spectral_info.v_primary = v_eig;
            spectral_info.w_generalized = w_gen;
        end

    else
        % =================================================================
        % Case 3: Distinct Real Eigenvalues
        % =================================================================
        l1 = (tr_A + sqrt(disc)) / 2.0;
        l2 = (tr_A - sqrt(disc)) / 2.0;
        [V_raw, ~] = eig(A_mat);
        v1 = V_raw(:, 1); if abs(v1(1)) > 1e-9, v1 = v1 / v1(1); end
        v2 = V_raw(:, 2); if abs(v2(1)) > 1e-9, v2 = v2 / v2(1); end
        
        u1_handle = @(t) exp(l1 * t) .* v1;
        u2_handle = @(t) exp(l2 * t) .* v2;
        
        spectral_info.case_type = "Distinct Real Eigenvalues (Node / Saddle)";
        spectral_info.eigenvalues = [l1; l2];
        spectral_info.v_primary = v1;
        spectral_info.w_generalized = v2;
    end

    %% 3. Match Initial Condition u(0) = a * u1(0) + b * u2(0)
    Basis_at_0 = [u1_handle(0), u2_handle(0)];
    ab_weights = Basis_at_0 \ u0;
    a_val = ab_weights(1);
    b_val = ab_weights(2);

    u_total_handle = @(t) a_val * u1_handle(t) + b_val * u2_handle(t);

    spectral_info.u1_handle = u1_handle;
    spectral_info.u2_handle = u2_handle;
    spectral_info.u_total_handle = u_total_handle;
    spectral_info.ab_weights = ab_weights;

    %% 4. Numerical Simulation via Control System Toolbox 'initial'
    t_eval = linspace(t_span(1), t_span(2), 600);
    [u_sim, ~] = initial(sys_ss, u0, t_eval);
    u_analytic = u_total_handle(t_eval).';

    max_err = max(max(abs(u_sim - u_analytic)));

    sim_data.t = t_eval;
    sim_data.u_sim = u_sim;
    sim_data.u_analytic = u_analytic;
    sim_data.verification_max_error = max_err;

    %% 5. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Complex & Repeated Eigen Solver: MIT 18.03 Lec 34\n');
    fprintf('====================================================\n');
    fprintf('Matrix A                : [ %+.2f, %+.2f ; %+.2f, %+.2f ]\n', ...
            A_mat(1,1), A_mat(1,2), A_mat(2,1), A_mat(2,2));
    fprintf('Invariants              : tr(A) = %+.4f, det(A) = %+.4f\n', tr_A, det_A);
    fprintf('Spectral Classification : %s\n', spectral_info.case_type);
    fprintf('Eigenvalues             : %+.4f %+.4fi,  %+.4f %+.4fi\n', ...
            real(spectral_info.eigenvalues(1)), imag(spectral_info.eigenvalues(1)), ...
            real(spectral_info.eigenvalues(2)), imag(spectral_info.eigenvalues(2)));
    fprintf('Basis Initial Vectors   : u1(0) = [%+.2f; %+.2f], u2(0) = [%+.2f; %+.2f]\n', ...
            Basis_at_0(1,1), Basis_at_0(2,1), Basis_at_0(1,2), Basis_at_0(2,2));
    fprintf('Superposition Weights   : a = %+.4f, b = %+.4f\n', a_val, b_val);
    fprintf('Verification Max Error  : %e (|ss initial - analytic|)\n', max_err);
    fprintf('====================================================\n\n');

    %% 6. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [100, 100, 1120, 460]);

        % Left Subplot: Phase Portrait
        subplot(1, 2, 1);
        hold on; grid on; box on;
        xline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        max_b = max([max(max(abs(u_sim))) * 1.2, 1.5]);
        [Xg, Yg] = meshgrid(linspace(-max_b, max_b, 15));
        quiver(Xg, Yg, A_mat(1,1)*Xg + A_mat(1,2)*Yg, A_mat(2,1)*Xg + A_mat(2,2)*Yg, ...
               'Color', [0.78 0.78 0.78], 'HandleVisibility', 'off');

        if isreal(spectral_info.v_primary)
            s_ray = linspace(-max_b, max_b, 100);
            vp = spectral_info.v_primary;
            plot(s_ray * vp(1), s_ray * vp(2), 'r--', 'LineWidth', 1.8, ...
                 'DisplayName', sprintf('Eigenline v = [%.1f; %.1f]', vp(1), vp(2)));
        end

        plot(u_sim(:, 1), u_sim(:, 2), 'b-', 'LineWidth', 2.2, 'DisplayName', 'Trajectory u(t)');
        plot(u0(1), u0(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Initial State u(0)');
        plot(u_sim(end, 1), u_sim(end, 2), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'DisplayName', 'Final State');

        xlabel('State x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State y(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Phase Portrait: %s', spectral_info.case_type), 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([-max_b, max_b]); ylim([-max_b, max_b]);

        % Right Subplot: Time-Domain Trajectories
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, u_sim(:, 1), 'b-', 'LineWidth', 2.0, 'DisplayName', 'State x(t) (ss initial)');
        plot(t_eval, u_sim(:, 2), 'r-', 'LineWidth', 2.0, 'DisplayName', 'State y(t) (ss initial)');
        plot(t_eval(1:30:end), u_analytic(1:30:end, 1), 'ko', 'MarkerSize', 5, 'DisplayName', 'Analytic x(t)');
        plot(t_eval(1:30:end), u_analytic(1:30:end, 2), 'go', 'MarkerSize', 5, 'DisplayName', 'Analytic y(t)');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title('Time-Domain Response & Analytical Verification', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end