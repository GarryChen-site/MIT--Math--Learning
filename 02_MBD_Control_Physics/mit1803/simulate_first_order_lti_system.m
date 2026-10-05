function [sim_results, sys_props] = simulate_first_order_lti_system(tau, gain, x0, u_handle, t_span)
% SIMULATE_FIRST_ORDER_LTI_SYSTEM Simulates an arbitrary first-order LTI system.
% Governing differential equation: tau * dx/dt + x = gain * u(t)
%
% % MBD Concept Mapping: Encapsulates standard first-order plant transfer function
%                      G(s) = gain / (tau*s + 1), performing time-domain response 
%                      decomposition into Zero-Input and Zero-State responses.
%
% Inputs (All optional with robust default values):
%   tau      - Time constant of the system in seconds (default: 0.001 s)
%   gain     - Steady-state static gain K (default: 1.0)
%   x0       - Initial state condition x(0) (default: 2.0)
%   u_handle - Function handle for input signal u(t) (default: Unit Step @ 5V)
%   t_span   - [t_start, t_end] simulation interval (default: [0, 0.005] s)
%
% Outputs:
%   sim_results - Structure containing arrays of [t, u, y_zir, y_zsr, y_total]
%   sys_props   - Structure containing poles, cutoff frequencies, and time constant

    arguments
        tau (1,1) double {mustBePositive} = 0.001
        gain (1,1) double {mustBeNumeric} = 1.0
        x0 (1,1) double {mustBeNumeric} = 2.0
        u_handle function_handle = @(t) 5.0 * ones(size(t))
        t_span (1,2) double {mustBeNumeric} = [0, 0.005]
    end

    %% 1. Analytical and Structural Properties Extraction
    sys_props.tau = tau;
    sys_props.gain = gain;
    sys_props.pole = -1 / tau;
    sys_props.omega_cutoff_rad = 1 / tau;
    sys_props.f_cutoff_hz = 1 / (2 * pi * tau);

    %% 2. Numerical Integration via High-Precision ODE Engine
    % Define the system state equation: dx/dt = (gain * u(t) - x) / tau
    ode_total = @(t, x) (gain * u_handle(t) - x) / tau;
    ode_zir   = @(t, x) (-x) / tau;
    ode_zsr   = @(t, x) (gain * u_handle(t) - x) / tau;

    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    
    % Uniform evaluation time grid
    t_eval = linspace(t_span(1), t_span(2), 1000);
    
    % Solve ZIR (x(0) = x0, u(t) = 0)
    [~, y_zir] = ode45(ode_zir, t_eval, x0, opts);
    
    % Solve ZSR (x(0) = 0, u(t) active)
    [~, y_zsr] = ode45(ode_zsr, t_eval, 0, opts);
    
    % Solve Total Response (x(0) = x0, u(t) active)
    [~, y_total] = ode45(ode_total, t_eval, x0, opts);

    % Package simulation output structure
    sim_results.t = t_eval;
    sim_results.u = u_handle(t_eval);
    sim_results.y_zir = y_zir';
    sim_results.y_zsr = y_zsr';
    sim_results.y_total = y_total';

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   First-Order LTI Dynamic Simulator (MBD Standard) \n');
    fprintf('====================================================\n');
    fprintf('System Equation         : %.4e * dx/dt + x = %.2f * u(t)\n', tau, gain);
    fprintf('Time Constant (tau)     : %.4e seconds\n', sys_props.tau);
    fprintf('Continuous Pole         : %.4e rad/s\n', sys_props.pole);
    fprintf('Cutoff Frequency        : %.2f Hz\n', sys_props.f_cutoff_hz);
    fprintf('Initial State Condition : x(0) = %.4f\n', x0);
    
    % Superposition error validation
    superposition_err = max(abs(sim_results.y_total - (sim_results.y_zir + sim_results.y_zsr)));
    fprintf('Superposition Error     : %e\n', superposition_err);
    fprintf('====================================================\n\n');

    %% 4. Standalone Demonstration Plot (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 850, 520]);
        hold on; grid on; box on;
        
        plot(sim_results.t, sim_results.u, 'k:', 'LineWidth', 1.5, ...
             'DisplayName', 'Input Signal u(t)');
        plot(sim_results.t, sim_results.y_zir, 'r--', 'LineWidth', 1.5, ...
             'DisplayName', 'Zero-Input Response (ZIR)');
        plot(sim_results.t, sim_results.y_zsr, 'b-.', 'LineWidth', 1.5, ...
             'DisplayName', 'Zero-State Response (ZSR)');
        plot(sim_results.t, sim_results.y_total, 'g-', 'LineWidth', 2.0, ...
             'DisplayName', 'Total Response (ZIR + ZSR)');

        xlabel('Time t (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('State Variable x(t)', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('LTI System Response: \\tau = %.2e s, Pole = %.1f rad/s', ...
              tau, sys_props.pole), 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 9);
        hold off;
    end
end