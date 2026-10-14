function [square_meta, sys_s_domain, sim_data] = solve_laplace_second_order_square_completion(plant_den, poly_input, x_init, t_span)
% SOLVE_LAPLACE_SECOND_ORDER_SQUARE_COMPLETION Solves 2nd-order LTI systems via Laplace & Completing the Square.
% Normalizes p(s) = a2*s^2 + a1*s + a0 into a2 * [ (s + sigma)^2 + kappa ] and resolves time response.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, impulse, pole, pzplot).
%
% % MBD Concept Mapping: Implements Completing the Square to extract modal damping decay sigma
% %                      and damped oscillation frequency wd directly from transfer function W(s) = 1/p(s).
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 28 Section [4] Benchmark):
%   plant_den  - Vector [a2, a1, a0] for a2*x'' + a1*x' + a0*x (default: [1.0, 2.0, 2.0])
%   poly_input - Polynomial forcing coefficients [b0, b1, ...] for b0 + b1*t + ...
%                If empty [], solves for Unit Impulse Response q(t) = delta(t) (default: [])
%   x_init     - Initial states [x(0+), x'(0+)] for regular ODEs (default: [0.0, 0.0])
%   t_span     - [t_start, t_end] simulation interval (default: [0, 6.0])
%
% Outputs:
%   square_meta  - Struct containing completed square parameters: sigma, kappa, wd, and regime
%   sys_s_domain - Struct containing W_tf (Transfer Function) and X_total_tf (Response Transform)
%   sim_data     - Struct containing time vector, impulse-inverted trajectory, and ode45 verification

    arguments
        plant_den (1,3) double {mustBeNumeric} = [1.0, 2.0, 2.0]
        poly_input double {mustBeNumeric} = []
        x_init (1,2) double {mustBeNumeric} = [0.0, 0.0]
        t_span (1,2) double {mustBeNumeric} = [0, 6.0]
    end

    %% 1. Completing the Square Algebraic Decomposition
    % p(s) = a2 * (s^2 + (a1/a2)*s + a0/a2) = a2 * [ (s + sigma)^2 + kappa ]
    a2 = plant_den(1); a1 = plant_den(2); a0 = plant_den(3);
    if abs(a2) < 1e-12
        error('Leading coefficient a2 must be non-zero for a second-order system.');
    end

    sigma = a1 / (2 * a2);
    kappa = (a0 / a2) - sigma^2;

    square_meta.a2 = a2;
    square_meta.sigma = sigma;
    square_meta.kappa = kappa;

    if kappa > 1e-12
        square_meta.wd = sqrt(kappa);
        square_meta.regime = "Underdamped: (s + sigma)^2 + wd^2 -> exp(-sigma*t)*sin(wd*t)";
        square_meta.formula_str = sprintf('%.2f * [ (s %+.2f)^2 + (%.4f)^2 ]', a2, sigma, square_meta.wd);
    elseif abs(kappa) <= 1e-12
        square_meta.wd = 0.0;
        square_meta.regime = "Critically Damped: (s + sigma)^2 -> t*exp(-sigma*t)";
        square_meta.formula_str = sprintf('%.2f * (s %+.2f)^2', a2, sigma);
    else
        square_meta.wd = sqrt(-kappa);
        square_meta.regime = "Overdamped: (s + sigma)^2 - alpha^2 -> exp(-sigma*t)*sinh(alpha*t)";
        square_meta.formula_str = sprintf('%.2f * [ (s %+.2f)^2 - (%.4f)^2 ]', a2, sigma, square_meta.wd);
    end

    %% 2. Transfer Function W(s) & Total Response X(s) Construction
    W_tf = tf(1.0, plant_den);
    sys_s_domain.W_tf = W_tf;

    x0 = x_init(1); v0 = x_init(2);
    ic_poly = [a2 * x0, a2 * v0 + a1 * x0]; % From L[a2*x''_r + a1*x'_r]

    if isempty(poly_input)
        % Case A: Unit Impulse Response w(t) under rest initial conditions (q(t) = delta(t))
        X_total_tf = W_tf;
        is_impulse_mode = true;
    else
        % Case B: Polynomial forcing q(t) = sum_{k=0}^N b_k * t^k with initial conditions
        is_impulse_mode = false;
        N_deg = length(poly_input) - 1;
        
        % Common denominator for Q(s): s^(N_deg + 1)
        % Q(s) = sum_{k=0}^N b_k * k! / s^(k+1) = [ sum_{k=0}^N b_k * k! * s^(N - k) ] / s^(N + 1)
        q_num_vec = zeros(1, N_deg + 1);
        for k = 0:N_deg
            q_num_vec(k + 1) = poly_input(k + 1) * factorial(k);
        end
        
        den_s_pow = [1.0, zeros(1, N_deg + 1)]; % s^(N + 1)
        den_common = conv(plant_den, den_s_pow);
        
        % Total numerator: q_num_vec + ic_poly * s^(N + 1)
        ic_shifted = conv(ic_poly, den_s_pow);
        max_len = max(length(q_num_vec), length(ic_shifted));
        q_pad  = [zeros(1, max_len - length(q_num_vec)), q_num_vec];
        ic_pad = [zeros(1, max_len - length(ic_shifted)), ic_shifted];
        
        num_common = q_pad + ic_pad;
        X_total_tf = tf(num_common, den_common);
    end

    sys_s_domain.X_total_tf = X_total_tf;

    %% 3. Time-Domain Inversion via 'impulse' & Cross-Check via ode45
    t_eval = linspace(t_span(1), t_span(2), 800);
    [x_laplace, ~] = impulse(X_total_tf, t_eval);

    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    if is_impulse_mode
        % Equivalent homogeneous IVP: x(0+) = 0, x'(0+) = 1/a2
        ode_fun = @(t, X) [X(2); (-a1 * X(2) - a0 * X(1)) / a2];
        [~, X_ode] = ode45(ode_fun, t_eval, [0.0; 1.0 / a2], opts);
    else
        q_handle = @(t) polyval(fliplr(poly_input), t);
        ode_fun = @(t, X) [X(2); (q_handle(t) - a1 * X(2) - a0 * X(1)) / a2];
        [~, X_ode] = ode45(ode_fun, t_eval, [x0; v0], opts);
    end
    x_ode = X_ode(:, 1);

    max_err = max(abs(x_laplace(:) - x_ode(:)));

    sim_data.t = t_eval;
    sim_data.x_laplace = x_laplace;
    sim_data.x_ode = x_ode;
    sim_data.verification_max_error = max_err;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Laplace 2nd-Order Solver: MIT 18.03 Lecture 28   \n');
    fprintf('====================================================\n');
    fprintf('Characteristic Poly p(s): %.2fs^2 %+.2fs %+.2f\n', a2, a1, a0);
    fprintf('Completed Square Form   : %s\n', square_meta.formula_str);
    fprintf('Damping Regime          : %s\n', square_meta.regime);
    fprintf('Transfer Function W(s)  :\n');
    display(W_tf);
    fprintf('Verification Max Error  : %e (|L^-1{X(s)} - ode45|)\n', max_err);
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 1050, 460]);

        % Left Subplot: Time-Domain Response & Exponential Envelope
        subplot(1, 2, 1);
        hold on; grid on; box on;
        plot(t_eval, x_laplace, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Laplace Solution L^{-1}\{X(s)\}');
        plot(t_eval(1:30:end), x_ode(1:30:end), 'ro', 'MarkerSize', 5, 'DisplayName', 'ode45 Equivalent Check');

        if is_impulse_mode && kappa > 0
            env_amp = 1.0 / (a2 * square_meta.wd);
            plot(t_eval, env_amp * exp(-sigma * t_eval), 'k--', 'LineWidth', 1.2, ...
                 'DisplayName', sprintf('Envelope \\pm %.2fe^{-%.2ft}', env_amp, sigma));
            plot(t_eval, -env_amp * exp(-sigma * t_eval), 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        end

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Response Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Time Response: %s', square_meta.formula_str), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);

        % Right Subplot: Pole-Zero Geometry in the s-Plane
        subplot(1, 2, 2);
        hold on; grid on; box on;
        sys_poles = pole(W_tf);
        xline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        yline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        xline(-sigma, 'r--', 'LineWidth', 1.4, ...
              'DisplayName', sprintf('Completed Square Axis: Re(s) = -\\sigma = %+.2f', -sigma));
        plot(real(sys_poles), imag(sys_poles), 'bx', 'MarkerSize', 10, 'LineWidth', 2.2, ...
             'DisplayName', 'System Poles (Roots of p(s))');

        xlabel('Real Axis Re(s)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Imaginary Axis Im(s)', 'FontSize', 10, 'FontWeight', 'bold');
        title('s-Plane Geometry: Completing the Square', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([min(real(sys_poles)) - 1.5, max(0.5, max(real(sys_poles)) + 1.0)]);
        ylim([-max(abs(imag(sys_poles))) - 1.5, max(abs(imag(sys_poles))) + 1.5]);
        hold off;
    end
end