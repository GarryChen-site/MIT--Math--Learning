function [m_avg, step_details, sol_trajectory] = solve_rk2_heun_step(f_ode, x0, y0, h, num_steps)
% SOLVE_RK2_HEUN_STEP Performs step-by-step RK2 (Heun's Method) integration.
%
% % Mathematical Tool Application: Second-order Runge-Kutta predictor-corrector.
%                                Demonstrates stage evaluations and slope averaging.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 2 Question 2.1):
%   f_ode     - Function handle @(x, y) defining dy/dx = f(x, y)
%   x0        - Initial independent variable value (default: 0)
%   y0        - Initial dependent state value (default: -1)
%   h         - Integration step size (default: 0.5)
%   num_steps - Number of integration steps to advance (default: 2)
%
% Outputs:
%   m_avg          - Effective average slope used in the first integration step
%   step_details   - Structure containing predictor/corrector points and slopes
%   sol_trajectory - Struct with array of evaluated state points [x, y]

    arguments
        f_ode function_handle = @(x, y) y^2 - x
        x0 (1,1) double {mustBeNumeric} = 0
        y0 (1,1) double {mustBeNumeric} = -1
        h (1,1) double {mustBePositive} = 0.5
        num_steps (1,1) double {mustBeInteger, mustBePositive} = 2
    end

    %% 1. Step 1 Detailed Breakdown (Solves MIT 18.03 Q2.1)
    % Initial slope polling (Predictor stage)
    m0 = f_ode(x0, y0);
    
    % Strut endpoint computation (Euler forward projection)
    x_predict = x0 + h;
    y_predict = y0 + h * m0;
    
    % Secondary slope polling (Corrector stage)
    m1 = f_ode(x_predict, y_predict);
    
    % Heun averaged slope
    m_avg = 0.5 * (m0 + m1);
    
    % Step 1 finalized state
    x1 = x_predict;
    y1 = y0 + h * m_avg;
    
    % Package step 1 diagnostic structure
    step_details.x0 = x0;
    step_details.y0 = y0;
    step_details.m0 = m0;
    step_details.x_predict = x_predict;
    step_details.y_predict = y_predict;
    step_details.m1 = m1;
    step_details.m_avg = m_avg;
    step_details.x1 = x1;
    step_details.y1 = y1;

    %% 2. Full Multi-Step Propagation
    x_arr = zeros(1, num_steps + 1);
    y_arr = zeros(1, num_steps + 1);
    x_arr(1) = x0;
    y_arr(1) = y0;
    
    for k = 1:num_steps
        curr_x = x_arr(k);
        curr_y = y_arr(k);
        
        slope_start = f_ode(curr_x, curr_y);
        y_temp = curr_y + h * slope_start;
        slope_end = f_ode(curr_x + h, y_temp);
        
        avg_slope = 0.5 * (slope_start + slope_end);
        
        x_arr(k+1) = curr_x + h;
        y_arr(k+1) = curr_y + h * avg_slope;
    end
    
    sol_trajectory.x = x_arr;
    sol_trajectory.y = y_arr;

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   MIT 18.03 Lec 2: RK2 (Heun) Step Evaluation     \n');
    fprintf('====================================================\n');
    fprintf('Initial State           : (x0, y0) = (%.4f, %.4f)\n', x0, y0);
    fprintf('Step Size (h)           : %.4f\n', h);
    fprintf('Predictor Slope (m0)    : %.4f\n', m0);
    fprintf('Euler Strut Point       : (x_pred, y_pred) = (%.4f, %.4f)\n', x_predict, y_predict);
    fprintf('Corrector Slope (m1)    : %.4f\n', m1);
    fprintf('----------------------------------------------------\n');
    fprintf('Heun Average Slope (m)  : %.4f (Analytical: 3/8 = 0.375)\n', m_avg);
    fprintf('Step 1 Finalized State  : y(%.4f) = %.4f\n', x1, y1);
    fprintf('====================================================\n\n');

    %% 4. Demonstration Plot (Active on direct F5 invocation)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [200, 200, 800, 500]);
        hold on; grid on; box on;
        
        % Plot Euler strut (predictor line)
        plot([x0, x_predict], [y0, y_predict], 'r--', 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Euler Strut (Slope m0 = %.2f)', m0));
        
        % Plot RK2 finalized step segment
        plot([x0, x1], [y0, y1], 'b-', 'LineWidth', 2.0, ...
             'DisplayName', sprintf('RK2 Heun Step (Slope m_{avg} = %.3f)', m_avg));
         
        % Plot full multi-step trajectory
        plot(sol_trajectory.x, sol_trajectory.y, 'ko-', 'LineWidth', 1.2, ...
             'DisplayName', 'RK2 Multi-Step Trajectory');

        % Highlight key stage points
        plot(x0, y0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'Initial State (0, -1)');
        plot(x_predict, y_predict, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'DisplayName', 'Predictor Point (0.5, -0.5)');
        plot(x1, y1, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'Corrected State (0.5, -0.8125)');

        xlabel('Independent Variable x', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('State Variable y', 'FontSize', 11, 'FontWeight', 'bold');
        title('Geometric Construction of RK2 / Heun Step (MIT 18.03 Q2.1)', 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 9);
        hold off;
    end
end