function [sol_analytic, sim_num] = solve_second_order_homogeneous(m, b, k, x0, v0, t_span)
% SOLVE_SECOND_ORDER_HOMOGENEOUS Solves m*x'' + b*x' + k*x = 0 with initial conditions.
%
% % Required Toolboxes:
% %   - Base MATLAB (Zero external dependency).
%
% % MBD Concept Mapping: Demonstrates characteristic polynomial root extraction (System Poles),
% %                      modal superposition, and state-space trajectory evolution [x; x_dot].
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 11 Benchmark):
%   m      - Mass / Inertia coefficient (default: 1.0)
%   b      - Damping coefficient (default: 5.0)
%   k      - Spring / Stiffness coefficient (default: 4.0)
%   x0     - Initial position state x(0) (default: 2.0)
%   v0     - Initial velocity state x'(0) (default: -5.0)
%   t_span - [t_start, t_end] evaluation time interval (default: [0, 4.0])
%
% Outputs:
%   sol_analytic - Struct containing roots, coefficients, and analytical function handles:
%                  .x_handle, .v_handle, .poles, .damping_type, .c1, .c2
%   sim_num      - Struct containing numerical ode45 verification data [t, x, v]

    arguments
        m (1,1) double {mustBePositive} = 1.0
        b (1,1) double {mustBeNonnegative} = 5.0
        k (1,1) double {mustBeNonnegative} = 4.0
        x0 (1,1) double {mustBeNumeric} = 2.0
        v0 (1,1) double {mustBeNumeric} = -5.0
        t_span (1,2) double {mustBeNumeric} = [0, 4.0]
    end

    %% 1. Characteristic Polynomial Roots & Damping Classification
    % Characteristic equation: m*s^2 + b*s + k = 0
    discriminant = b^2 - 4 * m * k;
    sol_analytic.discriminant = discriminant;

    if discriminant > 1e-12
        % Case 1: Overdamped (Distinct Real Roots)
        sol_analytic.damping_type = "Overdamped (Distinct Real Roots)";
        r1 = (-b + sqrt(discriminant)) / (2 * m);
        r2 = (-b - sqrt(discriminant)) / (2 * m);
        sol_analytic.poles = [r1, r2];

        % Linear system: [1, 1; r1, r2] * [c1; c2] = [x0; v0]
        c2 = (v0 - r1 * x0) / (r2 - r1);
        c1 = x0 - c2;
        sol_analytic.c1 = c1;
        sol_analytic.c2 = c2;

        sol_analytic.x_handle = @(t) c1 * exp(r1 * t) + c2 * exp(r2 * t);
        sol_analytic.v_handle = @(t) c1 * r1 * exp(r1 * t) + c2 * r2 * exp(r2 * t);

    elseif abs(discriminant) <= 1e-12
        % Case 2: Critically Damped (Repeated Real Root)
        sol_analytic.damping_type = "Critically Damped (Repeated Root)";
        r = -b / (2 * m);
        sol_analytic.poles = [r, r];

        % Form: x(t) = (c1 + c2*t)*exp(r*t) -> x(0)=c1, x'(0)=r*c1 + c2
        c1 = x0;
        c2 = v0 - r * c1;
        sol_analytic.c1 = c1;
        sol_analytic.c2 = c2;

        sol_analytic.x_handle = @(t) (c1 + c2 * t) .* exp(r * t);
        sol_analytic.v_handle = @(t) (c2 + r * (c1 + c2 * t)) .* exp(r * t);

    else
        % Case 3: Underdamped / Undamped (Complex Conjugate Roots)
        if b == 0
            sol_analytic.damping_type = "Undamped Harmonic Oscillator";
        else
            sol_analytic.damping_type = "Underdamped (Complex Conjugate Roots)";
        end
        alpha = -b / (2 * m);
        omega_d = sqrt(4 * m * k - b^2) / (2 * m);
        sol_analytic.poles = [alpha + 1i * omega_d, alpha - 1i * omega_d];

        % Form: x(t) = exp(alpha*t) * (c1*cos(omega_d*t) + c2*sin(omega_d*t))
        % x(0) = c1 -> c1 = x0; x'(0) = alpha*c1 + omega_d*c2 -> c2 = (v0 - alpha*c1)/omega_d
        c1 = x0;
        c2 = (v0 - alpha * c1) / omega_d;
        sol_analytic.c1 = c1;
        sol_analytic.c2 = c2;

        sol_analytic.x_handle = @(t) exp(alpha * t) .* (c1 * cos(omega_d * t) + c2 * sin(omega_d * t));
        sol_analytic.v_handle = @(t) exp(alpha * t) .* ((alpha * c1 + omega_d * c2) * cos(omega_d * t) + ...
                                                        (alpha * c2 - omega_d * c1) * sin(omega_d * t));
    end

    %% 2. Numerical Integration via State-Space ODE Representation
    % State variables: X = [x; v] -> dX/dt = [v; (-b*v - k*x)/m]
    ode_state = @(t, X) [X(2); (-b * X(2) - k * X(1)) / m];
    t_eval = linspace(t_span(1), t_span(2), 600);
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    
    [t_num, X_num] = ode45(ode_state, t_eval, [x0; v0], opts);
    sim_num.t = t_num;
    sim_num.x = X_num(:, 1);
    sim_num.v = X_num(:, 2);

    % Verification error assessment
    x_analytic_eval = sol_analytic.x_handle(sim_num.t);
    max_error = max(abs(x_analytic_eval - sim_num.x));
    sol_analytic.verification_max_error = max_error;

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Second-Order ODE Solver: MIT 18.03 Lecture 11   \n');
    fprintf('====================================================\n');
    fprintf('Equation Parameters     : m = %.2f, b = %.2f, k = %.2f\n', m, b, k);
    fprintf('System Classification   : %s\n', sol_analytic.damping_type);
    fprintf('Characteristic Poles    : %+.4f %+.4fi,  %+.4f %+.4fi\n', ...
            real(sol_analytic.poles(1)), imag(sol_analytic.poles(1)), ...
            real(sol_analytic.poles(2)), imag(sol_analytic.poles(2)));
    fprintf('Initial State           : x(0) = %.4f, v(0) = %.4f\n', x0, v0);
    fprintf('Superposition Constants : c1 = %.4f, c2 = %.4f\n', sol_analytic.c1, sol_analytic.c2);
    fprintf('Numeric Check Max Error : %e\n', max_error);
    fprintf('====================================================\n\n');

    %% 4. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [120, 120, 1100, 480]);

        % Left Subplot: Time-Domain State Evolution
        subplot(1, 2, 1);
        hold on; grid on; box on;
        plot(t_eval, sol_analytic.x_handle(t_eval), 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Analytic Position x(t)');
        plot(t_eval, sol_analytic.v_handle(t_eval), 'r--', 'LineWidth', 1.5, ...
             'DisplayName', 'Analytic Velocity \dot{x}(t)');
        plot(sim_num.t(1:25:end), sim_num.x(1:25:end), 'ko', 'MarkerSize', 5, ...
             'DisplayName', 'ode45 Position Check');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Time Domain Response (%s)', sol_analytic.damping_type), ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);

        % Right Subplot: Phase Plane Trajectory [x vs x_dot]
        subplot(1, 2, 2);
        hold on; grid on; box on;
        xline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        plot(sim_num.x, sim_num.v, 'b-', 'LineWidth', 1.8, ...
             'DisplayName', 'Phase Trajectory (x, \dot{x})');
        plot(x0, v0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, ...
             'DisplayName', 'Initial State (x_0, v_0)');
        plot(sim_num.x(end), sim_num.v(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 7, ...
             'DisplayName', 'Final State');

        xlabel('Position State x', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Velocity State \dot{x}', 'FontSize', 10, 'FontWeight', 'bold');
        title('Phase Plane State Trajectory', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        hold off;
    end
end