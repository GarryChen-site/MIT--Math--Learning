function [sol_analytic, sim_num] = solve_newton_cooling_ramp(k, T0, R_ramp, x0, t_span)
% SOLVE_NEWTON_COOLING_RAMP Solves first-order Newton's cooling ODE under ramp ambient.
% Governing differential equation: dx/dt + k*x = k*(T0 + R_ramp*t)
%
% % Required Toolboxes: None (Base MATLAB only).
%
% % MBD Concept Mapping: First-order thermal system dynamics. Deconstructs total
%                      solution into steady-state forced response (particular solution)
%                      and decaying transient (homogeneous solution).
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 4 benchmark):
%   k       - Thermal coupling constant (1/s or 1/h), default: 1/3
%   T0      - Ambient initial temperature at t = 0, default: 60
%   R_ramp  - Ambient ramp heating rate (deg/time), default: 6
%   x0      - Initial core temperature x(0), default: 32
%   t_span  - [t_start, t_end] evaluation time interval, default: [0, 15]
%
% Outputs:
%   sol_analytic - Structure containing analytical function handles and metrics:
%                  .x_total_handle: Full analytical solution x(t)
%                  .x_p_handle    : Particular steady-state solution x_p(t)
%                  .x_h_handle    : Decaying transient homogeneous solution x_h(t)
%                  .C_constant    : Integration constant from initial condition
%                  .lag           : Steady-state lag behind ambient (deg)
%                  .tau           : System time constant (1/k)
%   sim_num      - Numerical cross-check trajectory via ode45 [t, x]

    arguments
        k (1,1) double {mustBePositive} = 1/3
        T0 (1,1) double {mustBeNumeric} = 60
        R_ramp (1,1) double {mustBeNumeric} = 6
        x0 (1,1) double {mustBeNumeric} = 32
        t_span (1,2) double {mustBeNumeric} = [0, 15]
    end

    %% 1. Analytical Derivation via Integrating Factor
    % ODE: dx/dt + k*x = k*T0 + k*R_ramp*t
    % Integrating factor u(t) = exp(k*t)
    % Particular solution (method of undetermined coefficients or integration by parts):
    % Assume x_p(t) = A*t + B
    % x_p' + k*x_p = A + k*(A*t + B) = (k*A)*t + (A + k*B) = (k*R_ramp)*t + k*T0
    % -> A = R_ramp
    % -> A + k*B = k*T0 -> B = T0 - R_ramp/k
    A_coeff = R_ramp;
    B_coeff = T0 - (R_ramp / k);
    
    sol_analytic.tau = 1 / k;
    sol_analytic.lag = R_ramp / k;
    sol_analytic.x_p_handle = @(t) B_coeff + A_coeff * t;

    % Transient homogeneous term: x_h(t) = exp(-k*t)
    % Initial condition enforcement: x(0) = B_coeff + C = x0 -> C = x0 - B_coeff
    C_const = x0 - B_coeff;
    sol_analytic.C_constant = C_const;
    sol_analytic.x_h_handle = @(t) C_const * exp(-k * t);
    sol_analytic.x_total_handle = @(t) sol_analytic.x_p_handle(t) + sol_analytic.x_h_handle(t);

    %% 2. High-Precision Numerical Verification via ode45
    t_eval = linspace(t_span(1), t_span(2), 600);
    T_ext_eval = @(t) T0 + R_ramp * t;
    ode_fun = @(t, x) k * (T_ext_eval(t) - x);
    
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    [t_num, x_num] = ode45(ode_fun, t_eval, x0, opts);
    
    sim_num.t = t_num;
    sim_num.x = x_num;

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Newton Cooling Solver (MIT 18.03 Lec 4 Benchmark)\n');
    fprintf('====================================================\n');
    fprintf('Coupling Constant (k)   : %.4f\n', k);
    fprintf('Time Constant (tau)     : %.4f hours\n', sol_analytic.tau);
    fprintf('Initial State Condition : x(0) = %.4f\n', x0);
    fprintf('Steady-State Part x_p(t): %.4f + %.4f * t\n', B_coeff, A_coeff);
    fprintf('Transient Part x_h(t)   : (%.4f) * exp(-%.4f * t)\n', C_const, k);
    fprintf('Velocity Lag (e_ss)     : %.4f deg\n', sol_analytic.lag);
    
    % Verification error calculation
    x_analytic_eval = sol_analytic.x_total_handle(sim_num.t);
    max_error = max(abs(x_analytic_eval - sim_num.x));
    fprintf('Max Discrepancy (Analytic vs ode45): %e\n', max_error);
    fprintf('====================================================\n\n');

    %% 4. Self-Executing Demo Plot (Active on F5)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [180, 180, 850, 520]);
        hold on; grid on; box on;
        
        plot(t_eval, T_ext_eval(t_eval), 'r--', 'LineWidth', 1.8, ...
             'DisplayName', 'Ambient Ramp T_{ext}(t)');
        plot(t_eval, sol_analytic.x_p_handle(t_eval), 'k-', 'LineWidth', 2.0, ...
             'DisplayName', 'Steady-State Response x_p(t)');
        plot(sim_num.t, sim_num.x, 'bo', 'MarkerIndices', 1:25:length(sim_num.t), ...
             'LineWidth', 1.2, 'DisplayName', 'Numeric ode45 Verification');
        plot(t_eval, sol_analytic.x_total_handle(t_eval), 'c-', 'LineWidth', 1.2, ...
             'DisplayName', 'Analytic Total Solution x(t)');

        xlabel('Time t', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('State Variable x(t)', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('Newton Cooling Response (k = %.3f, Lag = %.1f)', k, sol_analytic.lag), ...
              'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 9);
        hold off;
    end
end