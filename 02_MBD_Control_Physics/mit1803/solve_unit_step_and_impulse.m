function [sol_analytic, sys_tf, sim_data] = solve_unit_step_and_impulse(sys_order, coeffs, t_span)
% SOLVE_UNIT_STEP_AND_IMPULSE Solves unit step v(t) and unit impulse w(t) responses under rest ICs.
% Evaluates initial jump derivatives: v'(0+) and w'(0+) (Solves MIT 18.03 Q24.1, Q24.2, Q24.3).
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, step, impulse, pole).
%
% % MBD Concept Mapping: Demonstrates fundamental test inputs: Step v(t) for tracking metrics,
% %                      and Impulse w(t) as the transfer function kernel h(t) = L^-1{G(s)}.
% %                      Validates equivalent initial value conversion: m*Delta(v) = Impulse.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 24 2nd-Order Benchmark):
%   sys_order - Integer order of the differential equation: 1 or 2 (default: 2)
%   coeffs    - Vector of ODE coefficients:
%               If order = 1: [a1, a0] corresponding to a1*x' + a0*x = q(t)
%               If order = 2: [a2, a1, a0] corresponding to a2*x'' + a1*x' + a0*x = q(t)
%               (default: [1.0, 0.0, 4.0] -> m=1, b=0, k=4)
%   t_span    - [t_start, t_end] evaluation time interval (default: [0, 2*pi])
%
% Outputs:
%   sol_analytic - Struct containing analytical function handles and initial derivative values
%   sys_tf       - Continuous transfer function object G(s) from Control System Toolbox
%   sim_data     - Struct containing step and impulse response trajectory arrays

    arguments
        sys_order (1,1) double {mustBeMember(sys_order, [1, 2])} = 2
        coeffs double {mustBeNumeric} = [1.0, 0.0, 4.0]
        t_span (1,2) double {mustBeNumeric} = [0, 2 * pi]
    end

    %% 1. Transfer Function Construction via Control System Toolbox
    sys_tf = tf(1.0, coeffs);
    sol_analytic.poles = pole(sys_tf);
    sol_analytic.order = sys_order;

    %% 2. Analytical Closed-Form Evaluation & Jump Discontinuities
    if sys_order == 1
        % Equation: a1*x' + a0*x = q(t) -> x' + k*x = (1/a1)*q(t), k = a0/a1
        a1 = coeffs(1); a0 = coeffs(2);
        k = a0 / a1;
        
        % Step Response v(t) = (1/a0) * (1 - exp(-k*t))
        sol_analytic.v_handle = @(t) (1 / a0) * (1 - exp(-k * t));
        % Impulse Response w(t) = (1/a1) * exp(-k*t)
        sol_analytic.w_handle = @(t) (1 / a1) * exp(-k * t);
        
        % Initial derivatives at t = 0+
        sol_analytic.v_dot_0plus = 1 / a1;       % MIT 18.03 Q24.1: v'(0+) = 1/a1
        sol_analytic.w_dot_0plus = -a0 / (a1^2); % MIT 18.03 Q24.2: w'(0+) = -k/a1
        
    else
        % Equation: a2*x'' + a1*x' + a0*x = q(t)
        a2 = coeffs(1); a1 = coeffs(2); a0 = coeffs(3);
        
        disc = a1^2 - 4 * a2 * a0;
        
        if abs(a1) < 1e-12
            % Undamped case: m*x'' + k*x = q(t)
            wn = sqrt(a0 / a2);
            sol_analytic.v_handle = @(t) (1 / a0) * (1 - cos(wn * t));
            sol_analytic.w_handle = @(t) (1 / (a2 * wn)) * sin(wn * t);
        elseif disc < 0
            % Underdamped case
            alpha = -a1 / (2 * a2);
            wd = sqrt(4 * a2 * a0 - a1^2) / (2 * a2);
            phi = atan2(wd, -alpha);
            sol_analytic.v_handle = @(t) (1 / a0) * (1 - (sqrt(alpha^2 + wd^2)/wd)*exp(alpha*t).*sin(wd*t + phi));
            sol_analytic.w_handle = @(t) (1 / (a2 * wd)) * exp(alpha * t) .* sin(wd * t);
        else
            % Overdamped / Critically damped
            r = pole(sys_tf);
            sol_analytic.v_handle = @(t) (1/a0) + (r(2)*exp(r(1)*t) - r(1)*exp(r(2)*t))/(a0*(r(1)-r(2)));
            sol_analytic.w_handle = @(t) (exp(r(1)*t) - exp(r(2)*t))/(a2*(r(1)-r(2)));
        end
        
        % Initial derivatives at t = 0+
        sol_analytic.v_dot_0plus = 0.0;    % Step velocity is continuous at 0
        sol_analytic.w_dot_0plus = 1 / a2; % MIT 18.03 Q24.3: w'(0+) = 1/m
    end

    %% 3. Control System Toolbox Simulation (step & impulse)
    t_eval = linspace(t_span(1), t_span(2), 1000);
    [y_step, ~]    = step(sys_tf, t_eval);
    [y_impulse, ~] = impulse(sys_tf, t_eval);

    sim_data.t = t_eval;
    sim_data.y_step = y_step;
    sim_data.y_impulse = y_impulse;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Step & Impulse Solver: MIT 18.03 Lecture 24      \n');
    fprintf('====================================================\n');
    fprintf('System Order            : %d\n', sys_order);
    fprintf('Initial Tangent v''(0+)  : %.4f\n', sol_analytic.v_dot_0plus);
    fprintf('Initial Tangent w''(0+)  : %.4f\n', sol_analytic.w_dot_0plus);
    fprintf('LTI Transfer Function G(s):\n');
    display(sys_tf);
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 1050, 480]);

        % Left Subplot: Unit Step Response
        subplot(1, 2, 1);
        hold on; grid on; box on;
        plot(t_eval, y_step, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Step Response v(t)');
        t_tan = linspace(0, (t_span(2)-t_span(1))*0.15, 100);
        plot(t_tan, sol_analytic.v_dot_0plus * t_tan, 'k--', 'LineWidth', 1.2, ...
             'DisplayName', sprintf('Initial Slope v''(0+) = %.2f', sol_analytic.v_dot_0plus));
        
        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Unit Step Response v(t) (Order %d)', sys_order), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);

        % Right Subplot: Unit Impulse Response
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, y_impulse, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Impulse Response w(t)');
        plot(t_tan, y_impulse(1) + sol_analytic.w_dot_0plus * t_tan, 'k--', 'LineWidth', 1.2, ...
             'DisplayName', sprintf('Initial Slope w''(0+) = %.2f', sol_analytic.w_dot_0plus));

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Impulse Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Unit Impulse Response w(t) (Order %d)', sys_order), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end