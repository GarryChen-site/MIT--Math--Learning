function [sol_analytic, sys_metrics] = solve_second_order_vibrations(m, b, k, x0, v0, t_span)
% SOLVE_SECOND_ORDER_VIBRATIONS Analyzes unforced vibration system mx'' + bx' + kx = 0.
%
% % Required Toolboxes: 
% %   - Base MATLAB (Zero external dependency).
%
% % MBD Concept Mapping: Extracts continuous system damping ratio zeta, undamped natural
% %                      frequency w_n, damped frequency w_d, and pole locations.
% %                      Validates Real Solution Theorem by decomposing complex modes.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 12 Underdamped Benchmark):
%   m      - Mass coefficient (default: 1.0)
%   b      - Damping constant (default: 4.0)
%   k      - Spring constant (default: 5.0)
%   x0     - Initial position x(0) (default: 1.0)
%   v0     - Initial velocity x'(0) (default: 0.0)
%   t_span - [t_start, t_end] simulation interval (default: [0, 6.0])
%
% Outputs:
%   sol_analytic - Struct containing function handles:
%                  .x_handle, .v_handle, .envelope_handle
%   sys_metrics  - Struct containing damping ratio, frequencies, poles, and classification

    arguments
        m (1,1) double {mustBePositive} = 1.0
        b (1,1) double {mustBeNonnegative} = 4.0
        k (1,1) double {mustBeNonnegative} = 5.0
        x0 (1,1) double {mustBeNumeric} = 1.0
        v0 (1,1) double {mustBeNumeric} = 0.0
        t_span (1,2) double {mustBeNumeric} = [0, 6.0]
    end

    %% 1. Analytical System Characteristics
    omega_n = sqrt(k / m);
    zeta = b / (2 * sqrt(m * k));
    disc = (b / (2 * m))^2 - (k / m);
    
    sys_metrics.omega_n = omega_n;
    sys_metrics.zeta = zeta;
    sys_metrics.critical_damping_b = 2 * sqrt(m * k);

    if disc < -1e-12
        % Case 1: Underdamped (Complex Conjugate Roots)
        sys_metrics.regime = "Underdamped";
        alpha = -b / (2 * m);
        omega_d = sqrt(-disc);
        period_d = 2 * pi / omega_d;
        
        sys_metrics.poles = [alpha + 1i*omega_d, alpha - 1i*omega_d];
        sys_metrics.omega_d = omega_d;
        sys_metrics.period_d = period_d;
        sys_metrics.decay_time_constant = -1 / alpha;

        % Enforcing x(0) = c1 = x0; x'(0) = alpha*c1 + omega_d*c2 = v0
        c1 = x0;
        c2 = (v0 - alpha * c1) / omega_d;
        A_polar = sqrt(c1^2 + c2^2);
        phi_polar = atan2(c2, c1);

        sol_analytic.x_handle = @(t) exp(alpha * t) .* (c1 * cos(omega_d * t) + c2 * sin(omega_d * t));
        sol_analytic.v_handle = @(t) exp(alpha * t) .* ((alpha * c1 + omega_d * c2) * cos(omega_d * t) + ...
                                                        (alpha * c2 - omega_d * c1) * sin(omega_d * t));
        sol_analytic.envelope_handle = @(t) A_polar * exp(alpha * t);
        sol_analytic.c_coefficients = [c1, c2];
        sol_analytic.polar_A_phi = [A_polar, phi_polar];

    elseif abs(disc) <= 1e-12
        % Case 2: Critically Damped (Repeated Real Root)
        sys_metrics.regime = "Critically Damped";
        r = -b / (2 * m);
        sys_metrics.poles = [r, r];
        sys_metrics.omega_d = 0.0;
        sys_metrics.period_d = Inf;
        sys_metrics.decay_time_constant = -1 / r;

        c1 = x0;
        c2 = v0 - r * c1;

        sol_analytic.x_handle = @(t) (c1 + c2 * t) .* exp(r * t);
        sol_analytic.v_handle = @(t) (c2 + r * (c1 + c2 * t)) .* exp(r * t);
        sol_analytic.envelope_handle = @(t) abs(c1 + c2 * t) .* exp(r * t);
        sol_analytic.c_coefficients = [c1, c2];
        sol_analytic.polar_A_phi = [NaN, NaN];

    else
        % Case 3: Overdamped (Distinct Real Roots)
        sys_metrics.regime = "Overdamped";
        r1 = -b / (2 * m) + sqrt(disc);
        r2 = -b / (2 * m) - sqrt(disc);
        sys_metrics.poles = [r1, r2];
        sys_metrics.omega_d = 0.0;
        sys_metrics.period_d = Inf;
        sys_metrics.decay_time_constant = -1 / max(r1, r2);

        c2 = (v0 - r1 * x0) / (r2 - r1);
        c1 = x0 - c2;

        sol_analytic.x_handle = @(t) c1 * exp(r1 * t) + c2 * exp(r2 * t);
        sol_analytic.v_handle = @(t) c1 * r1 * exp(r1 * t) + c2 * r2 * exp(r2 * t);
        sol_analytic.envelope_handle = @(t) abs(c1 * exp(r1 * t) + c2 * exp(r2 * t));
        sol_analytic.c_coefficients = [c1, c2];
        sol_analytic.polar_A_phi = [NaN, NaN];
    end

    %% 2. Numerical Integration Verification via ode45
    t_eval = linspace(t_span(1), t_span(2), 600);
    ode_fun = @(t, X) [X(2); -(k/m)*X(1) - (b/m)*X(2)];
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    
    [t_num, X_num] = ode45(ode_fun, t_eval, [x0; v0], opts);
    x_analytic_eval = sol_analytic.x_handle(t_num);
    max_error = max(abs(x_analytic_eval - X_num(:, 1)));
    sys_metrics.verification_max_error = max_error;

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Second-Order Vibration Solver: MIT 18.03 Lec 12  \n');
    fprintf('====================================================\n');
    fprintf('Physical System         : %.2f*x'''' + %.2f*x'' + %.2f*x = 0\n', m, b, k);
    fprintf('Damping Classification  : %s (Damping Ratio zeta = %.4f)\n', sys_metrics.regime, zeta);
    fprintf('Natural Frequency (w_n) : %.4f rad/s\n', omega_n);
    fprintf('Damped Frequency (w_d)  : %.4f rad/s\n', sys_metrics.omega_d);
    fprintf('Oscillation Period (P_d): %.4f s\n', sys_metrics.period_d);
    fprintf('System Poles (s)        : %+.4f %+.4fi,  %+.4f %+.4fi\n', ...
            real(sys_metrics.poles(1)), imag(sys_metrics.poles(1)), ...
            real(sys_metrics.poles(2)), imag(sys_metrics.poles(2)));
    fprintf('Verification Max Error  : %e\n', max_error);
    fprintf('====================================================\n\n');

    %% 4. Standalone Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [120, 120, 1050, 480]);

        % Left Subplot: Time-Domain State & Envelope
        subplot(1, 2, 1);
        hold on; grid on; box on;
        
        plot(t_eval, sol_analytic.x_handle(t_eval), 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Analytic Position x(t)');
        plot(t_eval, sol_analytic.envelope_handle(t_eval), 'r--', 'LineWidth', 1.4, ...
             'DisplayName', 'Decay Envelope +$|A|e^{-\sigma t}$');
        if sys_metrics.regime == "Underdamped"
            plot(t_eval, -sol_analytic.envelope_handle(t_eval), 'r--', 'LineWidth', 1.4, ...
                 'HandleVisibility', 'off');
        end
        plot(t_num(1:25:end), X_num(1:25:end, 1), 'ko', 'MarkerSize', 5, ...
             'DisplayName', 'ode45 Numerical Check');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Position x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Time Domain: %s (\\zeta = %.2f, \\omega_d = %.2f)', ...
              sys_metrics.regime, zeta, sys_metrics.omega_d), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5, 'Interpreter','latex');
        xlim([t_span(1), t_span(2)]);

        % Right Subplot: Phase Plane Trajectory
        subplot(1, 2, 2);
        hold on; grid on; box on;
        xline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        plot(X_num(:, 1), X_num(:, 2), 'b-', 'LineWidth', 1.8, ...
             'DisplayName', 'State Trajectory (x, \dot{x})');
        plot(x0, v0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, ...
             'DisplayName', 'Initial State (x_0, v_0)');
        plot(0, 0, 'rx', 'MarkerSize', 8, 'LineWidth', 1.5, ...
             'DisplayName', 'Equilibrium Origin');

        xlabel('Position State x', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Velocity State $\dot{x}$', 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter','latex');
        title('Phase Plane Orbit (Inward Spiral)', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5, 'Interpreter','latex');
        hold off;
    end
end