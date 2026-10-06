function [analysis_res, sim_data] = solve_autonomous_bifurcation_harvesting(a, h, y0, t_span)
% SOLVE_AUTONOMOUS_BIFURCATION_HARVESTING Analyzes autonomous ODE dy/dt = (a - y)*y - h.
%
% % Required Toolboxes: 
% %   - Base MATLAB (Zero external dependency).
%
% % MBD Concept Mapping: Demonstrates operating point extraction, local Jacobian
% %                      linearization, saddle-node bifurcation boundaries, and
% %                      basin of attraction thresholds.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 8 Benchmark):
%   a      - Capacity scale parameter, default: 1.3
%   h      - Constant harvesting rate, default: 0.25
%   y0     - Initial population state y(0), default: 0.5
%   t_span - [t_start, t_end] simulation interval, default: [0, 12.0]
%
% Outputs:
%   analysis_res - Struct containing bifurcation status, equilibria, and poles
%   sim_data     - Struct containing simulated trajectory arrays [t, y]

    arguments
        a (1,1) double {mustBeNumeric} = 1.3
        h (1,1) double {mustBeNonnegative} = 0.25
        y0 (1,1) double {mustBeNumeric} = 0.5
        t_span (1,2) double {mustBeNumeric} = [0, 12.0]
    end

    %% 1. Bifurcation and Operating Point Analysis
    discriminant = a^2 - 4 * h;
    a_crit = 2 * sqrt(h);
    
    analysis_res.a_parameter = a;
    analysis_res.h_parameter = h;
    analysis_res.a_critical = a_crit;
    analysis_res.inflection_point = a / 2;

    if discriminant < -1e-12
        analysis_res.bifurcation_state = "Subcritical: No Equilibrium (Definite Collapse)";
        analysis_res.equilibria = [];
        analysis_res.poles = [];
        analysis_res.stability = [];
    elseif abs(discriminant) <= 1e-12
        analysis_res.bifurcation_state = "Critical: Single Semistable Point (Saddle-Node Point)";
        y_star = a / 2;
        analysis_res.equilibria = y_star;
        analysis_res.poles = 0.0;
        analysis_res.stability = "Semistable";
    else
        analysis_res.bifurcation_state = "Supercritical: Two Distinct Operating Points";
        y_high = (a + sqrt(discriminant)) / 2;
        y_low  = (a - sqrt(discriminant)) / 2;
        
        analysis_res.equilibria = [y_high, y_low];
        
        pole_high = a - 2 * y_high;
        pole_low  = a - 2 * y_low;
        analysis_res.poles = [pole_high, pole_low];
        analysis_res.stability = ["Asymptotically Stable (Attractor)", "Unstable (Repeller / Threshold)"];
    end

    %% 2. Numerical Integration with Safe Event Guard
    ode_fun = @(t, y) (a - y) * y - h;
    
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Events', @local_event_extinction);
    
    [t_sol, y_sol, te, ~, ~] = ode45(ode_fun, t_span, y0, opts);
    
    sim_data.t = t_sol;
    sim_data.y = y_sol;
    sim_data.crashed = ~isempty(te);

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Autonomous System Analysis: MIT 18.03 Lec 8      \n');
    fprintf('====================================================\n');
    fprintf('Model Formulation      : dy/dt = (%.2f - y)*y - %.2f\n', a, h);
    fprintf('Critical Parameter (a) : %.4f (Saddle-Node Threshold)\n', a_crit);
    fprintf('Bifurcation Regime     : %s\n', analysis_res.bifurcation_state);
    fprintf('Point of Inflection    : y = %.4f (Where g''(y) = 0)\n', analysis_res.inflection_point);
    
    if ~isempty(analysis_res.equilibria)
        for idx = 1:length(analysis_res.equilibria)
            fprintf('  Equilibrium y*_%d      : %8.4f | Pole: %+.4f (%s)\n', ...
                    idx, analysis_res.equilibria(idx), analysis_res.poles(idx), analysis_res.stability(idx));
        end
    end
    fprintf('Simulation from y(0)=%.2f: Final y(t_end) = %.4f (Crashed: %d)\n', ...
            y0, y_sol(end), sim_data.crashed);
    fprintf('====================================================\n\n');

    %% 4. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 950, 480]);

        subplot(1, 2, 1);
        hold on; grid on; box on;
        y_grid = linspace(-0.2, a + 0.5, 300);
        g_vals = (a - y_grid) .* y_grid - h;
        
        plot(y_grid, g_vals, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Phase Velocity g(y) = dy/dt');
        yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        
        if ~isempty(analysis_res.equilibria)
            plot(analysis_res.equilibria(1), 0, 'bo', 'MarkerFaceColor', 'b', ...
                 'MarkerSize', 8, 'DisplayName', 'Stable Node y^*');
            if length(analysis_res.equilibria) > 1
                plot(analysis_res.equilibria(2), 0, 'ro', 'MarkerFaceColor', 'w', ...
                     'MarkerSize', 8, 'DisplayName', 'Unstable Node y^*');
            end
        end
        
        xlabel('State Variable y', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Velocity dy/dt = g(y)', 'FontSize', 10, 'FontWeight', 'bold');
        title('Phase Line Geometry & Operating Points', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'south', 'FontSize', 8.5);

        subplot(1, 2, 2);
        hold on; grid on; box on;
        
        plot(sim_data.t, sim_data.y, 'k-', 'LineWidth', 2.0, ...
             'DisplayName', sprintf('State Trajectory y(t), y_0=%.2f', y0));
         
        if ~isempty(analysis_res.equilibria)
            yline(analysis_res.equilibria(1), 'b--', 'LineWidth', 1.2, ...
                  'DisplayName', 'Attractor Asymptote');
            if length(analysis_res.equilibria) > 1
                yline(analysis_res.equilibria(2), 'r:', 'LineWidth', 1.2, ...
                      'DisplayName', 'Repeller Boundary');
            end
        end

        xlabel('Time t (seconds / years)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Variable y(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title('Dynamic Time Response', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end

%% Local Event Function
function [value, isterminal, direction] = local_event_extinction(~, y)
    value = y;          % Trigger when y drops to 0
    isterminal = 1;     % Halt integration
    direction = -1;     % Falling direction
end