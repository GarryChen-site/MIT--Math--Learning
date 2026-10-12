function [sol_analytic, sys_tf, sim_data] = solve_step_impulse_dynamics(m, b, k, t_span)
% SOLVE_STEP_IMPULSE_DYNAMICS Evaluates Step & Impulse responses of 2nd-order LTI plants.
% Verifies the Generalized Derivative identity: d/dt(y_step) = y_impulse = h(t).
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, step, impulse, stepinfo, pole).
%
% % MBD Concept Mapping: Demonstrates standard test input evaluation in control engineering:
% %                      Step input u(t) for tracking/overshoot metrics,
% %                      and Delta input delta(t) for system impulse kernel identification h(t) = L^-1{G(s)}.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 23 Benchmark):
%   m      - Mass / Inertia coefficient (default: 1.0)
%   b      - Damping constant (default: 1.2)
%   k      - Spring constant (default: 9.0)
%   t_span - [t_start, t_end] evaluation time interval (default: [0, 6.0])
%
% Outputs:
%   sol_analytic - Struct containing damping ratio, natural frequency, and stepinfo metrics
%   sys_tf       - Continuous transfer function object G(s) from Control System Toolbox
%   sim_data     - Struct containing step, impulse, and numerical gradient arrays

    arguments
        m (1,1) double {mustBePositive} = 1.0
        b (1,1) double {mustBeNonnegative} = 1.2
        k (1,1) double {mustBePositive} = 9.0
        t_span (1,2) double {mustBeNumeric} = [0, 6.0]
    end

    %% 1. Plant Transfer Function via Control System Toolbox
    % Differential equation: m*x'' + b*x' + k*x = u(t)
    % Transfer function: G(s) = 1 / (m*s^2 + b*s + k)
    sys_tf = tf(1.0, [m, b, k]);
    
    wn = sqrt(k / m);
    zeta = b / (2 * sqrt(m * k));
    
    sol_analytic.omega_n = wn;
    sol_analytic.zeta = zeta;
    sol_analytic.poles = pole(sys_tf);

    %% 2. Time-Domain Evaluation via Control System Toolbox
    t_eval = linspace(t_span(1), t_span(2), 1200);
    [y_step, t_step]       = step(sys_tf, t_eval);
    [y_impulse, t_impulse] = impulse(sys_tf, t_eval);

    % Extract transient performance metrics
    info_metrics = stepinfo(sys_tf);
    sol_analytic.step_metrics = info_metrics;

    %% 3. Numerical Verification of Generalized Derivative
    dt = t_eval(2) - t_eval(1);
    dy_step_dt = gradient(y_step, dt);
    
    % Maximum error between differentiated step response and impulse response
    diff_error = max(abs(dy_step_dt - y_impulse));
    sol_analytic.derivative_max_error = diff_error;

    sim_data.t = t_eval;
    sim_data.y_step = y_step;
    sim_data.y_impulse = y_impulse;
    sim_data.dy_step_dt = dy_step_dt;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Step & Delta LTI Dynamics: MIT 18.03 Lec 23     \n');
    fprintf('====================================================\n');
    fprintf('Physical Model          : %.2f*x'''' + %.2f*x'' + %.2f*x = u(t)\n', m, b, k);
    fprintf('Natural Frequency (w_n) : %.4f rad/s\n', wn);
    fprintf('Damping Ratio (zeta)    : %.4f\n', zeta);
    fprintf('Step Steady-State Value : %.4f (Analytical: 1/k = %.4f)\n', ...
            y_step(end), 1/k);
    fprintf('Step Rise Time (tr)     : %.4f s\n', info_metrics.RiseTime);
    fprintf('Step Overshoot (Mp)     : %.2f%%\n', info_metrics.Overshoot);
    fprintf('Derivative Max Error    : %e (|d/dt(step) - impulse|)\n', diff_error);
    fprintf('Plant Model G(s)        :\n');
    display(sys_tf);
    fprintf('====================================================\n\n');

    %% 5. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 1050, 480]);

        % Left Subplot: Unit Step Response
        subplot(1, 2, 1);
        hold on; grid on; box on;
        plot(t_eval, y_step, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Step Response x_{step}(t)');
        yline(1/k, 'k--', 'LineWidth', 1.0, 'DisplayName', sprintf('Steady-State Limit (1/k = %.3f)', 1/k));
        plot(info_metrics.PeakTime, info_metrics.Peak, 'ro', 'MarkerFaceColor', 'r', ...
             'MarkerSize', 6, 'DisplayName', sprintf('Peak Overshoot (%.1f%%)', info_metrics.Overshoot));

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Displacement x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Unit Step Response (\\omega_n = %.1f, \\zeta = %.2f)', wn, zeta), ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'southeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);

        % Right Subplot: Impulse Response & Gradient Verification
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, y_impulse, 'r-', 'LineWidth', 1.8, ...
             'DisplayName', 'Impulse Response h(t) = L^{-1}\{G(s)\}');
        plot(t_eval(1:30:end), dy_step_dt(1:30:end), 'ko', 'MarkerSize', 5, ...
             'DisplayName', 'Numerical Derivative d/dt(x_{step})');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Impulse Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title('Verification: Generalized Derivative d/dt(Step) = Impulse', ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end