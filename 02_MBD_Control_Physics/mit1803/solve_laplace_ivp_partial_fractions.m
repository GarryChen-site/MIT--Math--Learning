function [sol_s_domain, modal_residues, sim_data] = solve_laplace_ivp_partial_fractions(plant_den, B_amp, r_exp, x_init, t_span)
% SOLVE_LAPLACE_IVP_PARTIAL_FRACTIONS Solves LTI IVP p(D)x = B*e^{r*t} via Laplace Transform.
% Performs Partial Fraction Expansion (Cover-up / Residue method) and separates ZSR & ZIR.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, impulse, pole).
%
% % MBD Concept Mapping: Demonstrates Transfer Function W(s) = 1/p(s), algebraic IVP resolution
% %                      via L[f'_r] = sF(s) - f(0+), and modal residue decomposition.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 27 Section [5] Benchmark):
%   plant_den - Polynomial coefficients of p(s): [a1, a0] or [a2, a1, a0] (default: [1.0, 3.0])
%   B_amp     - Exponential input amplitude B in B*exp(r*t) (default: 1.0)
%   r_exp     - Exponential input rate r in B*exp(r*t) (default: -1.0)
%   x_init    - Initial state at t = 0+: scalar [x0] for 1st-order, [x0, v0] for 2nd-order (default: 5.0)
%   t_span    - [t_start, t_end] simulation interval (default: [0, 4.0])
%
% Outputs:
%   sol_s_domain   - Struct containing tf objects: .W_transfer, .X_zsr, .X_zir, .X_total
%   modal_residues - Struct containing partial fraction residues and poles
%   sim_data       - Struct containing time trajectories and ode45 verification error

    arguments
        plant_den (1,:) double {mustBeNumeric} = [1.0, 3.0]
        B_amp (1,1) double {mustBeNumeric} = 1.0
        r_exp (1,1) double {mustBeNumeric} = -1.0
        x_init (1,:) double {mustBeNumeric} = 5.0
        t_span (1,2) double {mustBeNumeric} = [0, 4.0]
    end

    %% 1. Construct Transfer Function & Common-Denominator s-Domain Polynomials
    sys_order = length(plant_den) - 1;

    % System Transfer Function: W(s) = 1 / p(s)
    W_tf = tf(1.0, plant_den);

    % Construct initial condition polynomial IC_poly(s) from L[f'_r] = s*F(s) - f(0+)
    if sys_order == 1
        % a1*(s*X - x0) + a0*X = Q(s)  -->  p(s)*X = Q(s) + a1*x0
        a1 = plant_den(1);
        x0 = x_init(1);
        ic_poly_vec = a1 * x0; % Scalar polynomial (degree 0)
    elseif sys_order == 2
        % a2*(s^2*X - s*x0 - v0) + a1*(s*X - x0) + a0*X = Q(s)
        % p(s)*X = Q(s) + a2*x0*s + (a2*v0 + a1*x0)
        a2 = plant_den(1); a1 = plant_den(2);
        x0 = x_init(1);
        v0 = 0.0;
        if length(x_init) >= 2
            v0 = x_init(2);
        end
        ic_poly_vec = [a2 * x0, a2 * v0 + a1 * x0]; % Degree 1 polynomial in s
    else
        error('Supported system order is 1 or 2.');
    end

    % Exact Common Denominator for forced system: D_common(s) = p(s) * (s - r_exp)
    den_common = conv(plant_den, [1.0, -r_exp]);

    % Numerator for ZSR: N_zsr(s) = B_amp
    num_zsr = B_amp;

    % Numerator for Total Response: N_tot(s) = B_amp + IC_poly(s) * (s - r_exp)
    ic_times_sr = conv(ic_poly_vec, [1.0, -r_exp]);
    num_tot = ic_times_sr;
    num_tot(end) = num_tot(end) + B_amp; % Add constant B_amp to the lowest power of s

    % Create exact LTI Transfer Function objects without uncancelled pole-zero pairs
    X_zsr   = tf(num_zsr, den_common);
    X_zir   = tf(ic_poly_vec, plant_den);
    X_total = tf(num_tot, den_common);

    sol_s_domain.W_transfer = W_tf;
    sol_s_domain.X_zsr = X_zsr;
    sol_s_domain.X_zir = X_zir;
    sol_s_domain.X_total = X_total;

    %% 2. Partial Fraction Expansion (Guaranteed Matching Dimensions & Pole Order)
    [res_tot, poles_tot, k_direct] = residue(num_tot, den_common);
    [res_zsr, poles_zsr, ~]        = residue(num_zsr, den_common);

    modal_residues.total_residues = res_tot;
    modal_residues.total_poles = poles_tot;
    modal_residues.zsr_residues = res_zsr;
    modal_residues.zsr_poles = poles_zsr;
    modal_residues.direct_term = k_direct;

    %% 3. Time-Domain Inversion via Control System Toolbox 'impulse' & ode45 Check
    t_eval = linspace(t_span(1), t_span(2), 800);
    [x_zsr_t, ~]   = impulse(X_zsr, t_eval);
    [x_zir_t, ~]   = impulse(X_zir, t_eval);
    [x_total_t, ~] = impulse(X_total, t_eval);

    % Numerical cross-check via ode45
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    if sys_order == 1
        ode_fun = @(t, x) (B_amp * exp(r_exp * t) - plant_den(2) * x) / plant_den(1);
        [~, x_ode] = ode45(ode_fun, t_eval, x_init(1), opts);
    else
        ode_fun = @(t, X) [X(2); (B_amp * exp(r_exp * t) - plant_den(2)*X(2) - plant_den(3)*X(1)) / plant_den(1)];
        [~, X_ode] = ode45(ode_fun, t_eval, [x0; v0], opts);
        x_ode = X_ode(:, 1);
    end

    max_err = max(abs(x_total_t(:) - x_ode(:)));

    sim_data.t = t_eval;
    sim_data.x_zsr = x_zsr_t;
    sim_data.x_zir = x_zir_t;
    sim_data.x_total = x_total_t;
    sim_data.x_ode = x_ode;
    sim_data.verification_max_error = max_err;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Laplace IVP & Residue Solver: MIT 18.03 Lec 27   \n');
    fprintf('====================================================\n');
    fprintf('System Transfer Function W(s) = 1/p(s):\n');
    display(W_tf);
    fprintf('--- Partial Fractions of X_zsr(s) (Cover-up Method) ---\n');
    for idx = 1:length(poles_zsr)
        fprintf('  Pole s = %+.4f %+.4fi  -->  Residue = %+.4f %+.4fi\n', ...
                real(poles_zsr(idx)), imag(poles_zsr(idx)), real(res_zsr(idx)), imag(res_zsr(idx)));
    end
    fprintf('--- Partial Fractions of Total X(s) (ZSR + ZIR) ---\n');
    for idx = 1:length(poles_tot)
        fprintf('  Pole s = %+.4f %+.4fi  -->  Residue = %+.4f %+.4fi\n', ...
                real(poles_tot(idx)), imag(poles_tot(idx)), real(res_tot(idx)), imag(res_tot(idx)));
    end
    fprintf('Verification Max Error (|L^-1{X(s)} - ode45|) : %e\n', max_err);
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 1050, 460]);

        % Left Subplot: ZSR, ZIR, and Total Response Superposition
        subplot(1, 2, 1);
        hold on; grid on; box on;
        plot(t_eval, x_total_t, 'k-', 'LineWidth', 2.2, 'DisplayName', 'Total Response x(t)');
        plot(t_eval, x_zsr_t, 'b--', 'LineWidth', 1.6, 'DisplayName', 'Zero-State Response x_{zsr}(t)');
        plot(t_eval, x_zir_t, 'r-.', 'LineWidth', 1.6, 'DisplayName', 'Zero-Input Transient x_{zir}(t)');
        plot(t_eval(1:30:end), x_ode(1:30:end), 'go', 'MarkerSize', 5, 'DisplayName', 'ode45 Check');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Amplitude x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title('Time-Domain Inversion: x(t) = L^{-1}\{X_{zsr}(s) + X_{zir}(s)\}', ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);

        % Right Subplot: Modal Residue Bar Chart (Using categorical labels for robustness)
        subplot(1, 2, 2);
        hold on; grid on; box on;
        
        pole_labels = strings(length(poles_tot), 1);
        for k = 1:length(poles_tot)
            if abs(imag(poles_tot(k))) < 1e-6
                pole_labels(k) = sprintf('s = %.2f', real(poles_tot(k)));
            else
                pole_labels(k) = sprintf('s = %.1f%+.1fi', real(poles_tot(k)), imag(poles_tot(k)));
            end
        end
        
        bar_data = [real(res_zsr(:)), real(res_tot(:))];
        b_handle = bar(categorical(pole_labels, pole_labels), bar_data, 'grouped');
        b_handle(1).DisplayName = 'ZSR Residues (x(0^+)=0)';
        b_handle(2).DisplayName = 'Total Residues (with IC)';

        xlabel('System & Input Poles', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Real Part of Residue Weight R_k', 'FontSize', 10, 'FontWeight', 'bold');
        title('Partial Fraction Modal Weights (Cover-up Method)', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        hold off;
    end
end