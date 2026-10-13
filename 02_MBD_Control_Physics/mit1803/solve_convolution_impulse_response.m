function [sim_results, sys_tf] = solve_convolution_impulse_response(plant_input, q_func, t_span, dt_step)
% SOLVE_CONVOLUTION_IMPULSE_RESPONSE Evaluates LTI system response via Convolution Integral.
% Implements Duhamel's Superposition Integral: x(t) = \int_0^t w(t - tau) * q(tau) dtau.
% Cross-checks numerical convolution against Control System Toolbox 'lsim'.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, impulse, lsim, isproper).
%
% % MBD Concept Mapping: Demonstrates non-parametric kernel extraction w(t) = L^-1{G(s)},
% %                      the equivalence of convolution to continuous LTI state integration,
% %                      and impulse superposition for arbitrary external drives.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 25 Lake Phosphate Benchmark):
%   plant_input - Control System Toolbox 'tf' object OR denominator vector [a1, a0] (default: tf(1, [1, 1]))
%   q_func      - Function handle for input signal q(t) (default: @(t) 1.0 + cos(2.0*t))
%   t_span      - [t_start, t_end] evaluation horizon (default: [0, 8.0])
%   dt_step     - Discretization time step for Riemann sum (default: 0.005)
%
% Outputs:
%   sim_results - Struct containing time vector, kernel w(t), convolution output, lsim output, and error
%   sys_tf      - Continuous transfer function object G(s) from Control System Toolbox

    arguments
        plant_input = tf(1.0, [1.0, 1.0]) % G(s) = 1 / (s + 1)
        q_func = @(t) 1.0 + cos(2.0 * t)  % Lecture phosphate input q(t) = 1 + cos(2t)
        t_span (1,2) double {mustBeNumeric} = [0, 8.0]
        dt_step (1,1) double {mustBePositive} = 0.005
    end

    %% 1. Plant Transfer Function Standardization via Control System Toolbox
    if isa(plant_input, 'tf') || isa(plant_input, 'zpk') || isa(plant_input, 'ss')
        sys_tf = tf(plant_input);
    elseif isnumeric(plant_input)
        % Treated as polynomial coefficients of denominator: [a_n, ..., a_0]
        sys_tf = tf(1.0, plant_input);
    else
        error('Invalid plant input. Provide a Control System Toolbox LTI object or numeric vector.');
    end

    %% 2. Time Grid and Input Signal Generation
    t_eval = t_span(1):dt_step:t_span(2);
    q_eval = q_func(t_eval);

    %% 3. Unit Impulse Response Extraction via Control System Toolbox
    [w_kernel, ~] = impulse(sys_tf, t_eval);

    %% 4. Numerical Convolution Evaluation (Duhamel / Superposition Integral)
    % Emulating discrete Riemann sum approximation from lecture:
    % x(t) = sum_k w(t - k*dt) * q(k*dt) * dt
    conv_full = conv(q_eval, w_kernel) * dt_step;
    x_conv = conv_full(1:length(t_eval)); % Truncate to evaluation interval length

    %% 5. Verification via Control System Toolbox 'lsim'
    [x_lsim, ~] = lsim(sys_tf, q_eval, t_eval);

    % Maximum absolute error assessment
    max_err = max(abs(x_conv(:) - x_lsim(:)));

    sim_results.t = t_eval;
    sim_results.q_input = q_eval;
    sim_results.w_kernel = w_kernel;
    sim_results.x_convolution = x_conv;
    sim_results.x_lsim = x_lsim;
    sim_results.max_error = max_err;

    %% 6. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   LTI Convolution Solver: MIT 18.03 Lecture 25    \n');
    fprintf('====================================================\n');
    fprintf('Plant Transfer Function G(s):\n');
    display(sys_tf);
    fprintf('Simulation Horizon      : [%.2f, %.2f] seconds (dt = %.4f s)\n', ...
            t_span(1), t_span(2), dt_step);
    fprintf('Impulse Kernel at t=0+  : w(0+) = %.4f\n', w_kernel(1));
    fprintf('Verification Max Error  : %e (|x_{conv} - x_{lsim}|)\n', max_err);
    fprintf('====================================================\n\n');

    %% 7. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [120, 120, 1100, 480]);

        % Left Subplot: Kernel & Input Excitation
        subplot(1, 2, 1);
        hold on; grid on; box on;
        plot(t_eval, q_eval, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6, ...
             'DisplayName', 'Input Excitation q(t)');
        plot(t_eval, w_kernel, 'Color', [0.49 0.18 0.56], 'LineWidth', 1.8, ...
             'DisplayName', 'Impulse Kernel w(t) = L^{-1}\{G(s)\}');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Signal Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title('Convolution Operators: Input q(t) & Kernel w(t)', ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);

        % Right Subplot: Superposition Response Comparison
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, x_conv, 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Convolution Integral: (w * q)(t)');
        plot(t_eval(1:35:end), x_lsim(1:35:end), 'ro', 'MarkerSize', 5, ...
             'DisplayName', 'Control Toolbox lsim Check');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Response x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Dynamic Response Comparison (Max Error: %.2e)', max_err), ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'southeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end