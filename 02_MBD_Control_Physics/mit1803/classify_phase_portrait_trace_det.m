function [portrait_report, sys_ss] = classify_phase_portrait_trace_det(A_mat, t_horizon)
% CLASSIFY_PHASE_PORTRAIT_TRACE_DET Classifies 2x2 system u' = A*u on the Trace-Determinant Plane.
% Determines exact geometric portrait type (including marginal cases), spiral rotation direction,
% and Lyapunov stability category according to MIT 18.03 Lecture 35 rules.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (ss, initial, pole).
%
% % MBD Concept Mapping: Automates 2nd-order state-space stability & modal topology classification
% %                      using matrix invariants (tr A, det A) and structural inspection.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 35 Sec [3] Defective Sink):
%   A_mat     - 2x2 system matrix [a, b; c, d] (default: [-2, 1; -1, 0])
%   t_horizon - Simulation time horizon for phase plane portrait (default: 4.0)
%
% Outputs:
%   portrait_report - Struct with trace, det, discriminant, portrait_name, stability, and rotation
%   sys_ss          - State-Space object from Control System Toolbox

    arguments
        A_mat (2,2) double {mustBeNumeric} = [-2.0, 1.0; -1.0, 0.0]
        t_horizon (1,1) double {mustBePositive} = 4.0
    end

    %% 1. Compute Matrix Invariants & Eigenvalues via Control System Toolbox
    sys_ss = ss(A_mat, zeros(2, 1), eye(2), zeros(2, 1));
    
    tr_A  = trace(A_mat);
    det_A = det(A_mat);
    disc  = tr_A^2 - 4 * det_A; % tr(A)^2 - 4*det(A)
    tol   = 1e-10;

    eigenvals = pole(sys_ss);

    %% 2. Complete Phase Portrait & Stability Classification (Lecture 35 Rules)
    rotation_str = "N/A (Non-spiraling / Non-elliptical)";

    if det_A < -tol
        % Lower Half-Plane: det < 0 -> Real eigenvalues of opposite sign
        portrait_name = "Saddle Point";
        stability_cat = "Unstable";

    elseif abs(det_A) <= tol
        % Horizontal Axis: det = 0 -> At least one zero eigenvalue
        if tr_A < -tol
            portrait_name = "Comb (Line of Equilibria with Decaying Parallel Rays)";
            stability_cat = "Neutrally Stable";
        elseif tr_A > tol
            portrait_name = "Unstable Comb (Line of Equilibria with Growing Parallel Rays)";
            stability_cat = "Unstable";
        else
            if norm(A_mat, 'fro') < tol
                portrait_name = "Zero Matrix (Entire Plane is Static Equilibrium)";
                stability_cat = "Neutrally Stable";
            else
                portrait_name = "Parallel Shear Flow (Double Zero Eigenvalue, Defective)";
                stability_cat = "Unstable (Linear Secular Drift)";
            end
        end

    else
        % Upper Half-Plane: det > 0
        if disc < -tol
            % Above Critical Parabola: Complex Conjugate Eigenvalues
            if A_mat(2, 1) > tol
                rotation_str = "Counterclockwise (CCW, bottom-left entry c > 0)";
            else
                rotation_str = "Clockwise (CW, bottom-left entry c < 0)";
            end

            if abs(tr_A) <= tol
                portrait_name = "Center (Closed Elliptical Orbits)";
                stability_cat = "Neutrally Stable";
            elseif tr_A < -tol
                portrait_name = "Stable Spiral (Spiral Sink)";
                stability_cat = "Asymptotically Stable";
            else
                portrait_name = "Unstable Spiral (Spiral Source)";
                stability_cat = "Unstable";
            end

        elseif abs(disc) <= tol
            % On the Critical Parabola: det = tr^2 / 4 (Repeated Real Eigenvalue)
            lambda_rep = tr_A / 2.0;
            is_complete = norm(A_mat - lambda_rep * eye(2), 'fro') < tol;

            if is_complete
                if tr_A < 0
                    portrait_name = "Star Sink (Complete Repeated Eigenvalue)";
                    stability_cat = "Asymptotically Stable";
                else
                    portrait_name = "Star Source (Complete Repeated Eigenvalue)";
                    stability_cat = "Unstable";
                end
            else
                if tr_A < 0
                    portrait_name = "Defective Nodal Sink (Improper Node)";
                    stability_cat = "Asymptotically Stable";
                else
                    portrait_name = "Defective Nodal Source (Improper Node)";
                    stability_cat = "Unstable";
                end
            end

        else
            % Below Critical Parabola (and det > 0): Distinct Real Eigenvalues of Same Sign
            if tr_A < 0
                portrait_name = "Stable Node (Nodal Sink)";
                stability_cat = "Asymptotically Stable";
            else
                portrait_name = "Unstable Node (Nodal Source)";
                stability_cat = "Unstable";
            end
        end
    end

    portrait_report.A = A_mat;
    portrait_report.trace = tr_A;
    portrait_report.determinant = det_A;
    portrait_report.discriminant = disc;
    portrait_report.eigenvalues = eigenvals;
    portrait_report.portrait_name = portrait_name;
    portrait_report.stability_category = stability_cat;
    portrait_report.rotation_direction = rotation_str;

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Phase Portrait Classifier: MIT 18.03 Lecture 35  \n');
    fprintf('====================================================\n');
    fprintf('Matrix A                : [ %+.2f, %+.2f ; %+.2f, %+.2f ]\n', ...
            A_mat(1,1), A_mat(1,2), A_mat(2,1), A_mat(2,2));
    fprintf('Trace & Determinant     : tr(A) = %+.4f,  det(A) = %+.4f\n', tr_A, det_A);
    fprintf('Parabola Discriminant   : tr(A)^2 - 4*det(A) = %+.4f\n', disc);
    fprintf('Eigenvalues             : %+.4f %+.4fi,  %+.4f %+.4fi\n', ...
            real(eigenvals(1)), imag(eigenvals(1)), real(eigenvals(2)), imag(eigenvals(2)));
    fprintf('Phase Portrait Type     : %s\n', portrait_name);
    fprintf('Stability Category      : %s\n', stability_cat);
    fprintf('Rotation Direction      : %s\n', rotation_str);
    fprintf('====================================================\n\n');

    %% 4. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [100, 100, 1120, 460]);

        % Left Subplot: Location on the Trace-Determinant Plane
        subplot(1, 2, 1);
        hold on; grid on; box on;
        tr_lim = max(abs(tr_A) * 1.5, 3.0);
        det_lim = max(abs(det_A) * 1.5, (tr_lim^2)/4);

        patch([-tr_lim, 0, 0, -tr_lim], [0, 0, det_lim, det_lim], [0.85 0.96 0.85], ...
              'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', 'Asymptotically Stable Region');

        tr_grid = linspace(-tr_lim, tr_lim, 300);
        plot(tr_grid, (tr_grid.^2)/4, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Critical Parabola det = tr^2/4');
        xline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        yline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');

        plot(tr_A, det_A, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9, ...
             'DisplayName', sprintf('Current System (tr=%.2f, det=%.2f)', tr_A, det_A));

        xlabel('Trace tr(A)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Determinant det(A)', 'FontSize', 10, 'FontWeight', 'bold');
        title('System Location on Trace-Determinant Plane', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([-tr_lim, tr_lim]); ylim([-det_lim * 0.6, det_lim]);

        % Right Subplot: Simulated Phase Portrait Family via Control Toolbox 'initial'
        subplot(1, 2, 2);
        hold on; grid on; box on;
        xline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        [Xg, Yg] = meshgrid(linspace(-2, 2, 15));
        quiver(Xg, Yg, A_mat(1,1)*Xg + A_mat(1,2)*Yg, A_mat(2,1)*Xg + A_mat(2,2)*Yg, ...
               'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');

        % Plot real eigenlines if they exist
        [V_eig, D_eig] = eig(A_mat);
        if isreal(D_eig)
            s_line = linspace(-2.5, 2.5, 100);
            plot(s_line * V_eig(1,1), s_line * V_eig(2,1), 'r--', 'LineWidth', 1.6, ...
                 'DisplayName', sprintf('Eigenline 1 (\\lambda=%.2f)', D_eig(1,1)));
            if rank(V_eig) == 2 && abs(D_eig(1,1) - D_eig(2,2)) > tol
                plot(s_line * V_eig(1,2), s_line * V_eig(2,2), 'm--', 'LineWidth', 1.6, ...
                     'DisplayName', sprintf('Eigenline 2 (\\lambda=%.2f)', D_eig(2,2)));
            end
        end

        % Launch a ring of initial conditions around the origin
        theta_ring = linspace(0, 2*pi, 13);
        t_sim = linspace(0, t_horizon, 300);
        for ang = theta_ring(1:end-1)
            u0_k = 1.2 * [cos(ang); sin(ang)];
            [u_tr, ~] = initial(sys_ss, u0_k, t_sim);
            plot(u_tr(:, 1), u_tr(:, 2), 'b-', 'LineWidth', 1.4, 'HandleVisibility', 'off');
            plot(u0_k(1), u0_k(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
        plot(NaN, NaN, 'b-', 'LineWidth', 1.4, 'DisplayName', 'Trajectories u(t)');
        plot(NaN, NaN, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 5, 'DisplayName', 'Initial States u(0)');

        xlabel('State x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State y(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('%s (%s)', portrait_name, stability_cat), 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8);
        axis equal; xlim([-2, 2]); ylim([-2, 2]);
        hold off;
    end
end