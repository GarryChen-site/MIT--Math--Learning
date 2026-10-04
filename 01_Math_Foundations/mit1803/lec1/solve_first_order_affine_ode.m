function [y_analytic_func, sol_num] = solve_first_order_affine_ode(a, b, x0, y0, x_span)
% SOLVE_FIRST_ORDER_AFFINE_ODE Solves dy/dx = a*y + b analytically and numerically.
%
% Mathematical Tool Application: Solves first-order linear affine ODEs governing
%                                standard first-order dynamic responses.
%
% Inputs (All optional with defaults to MIT 18.03 Lec 1 Question 1):
%   a       - Linear state coefficient (System Pole = a), default: 2
%   b       - Exogenous constant forcing term, default: 1
%   x0      - Initial independent variable value, default: 0
%   y0      - Initial state condition y(x0) = y0, default: 0
%   x_span  - [x_start, x_end] integration interval, default: [0, 1.5]
%
% Outputs:
%   y_analytic_func - Function handle of the symbolic closed-form solution
%   sol_num         - Structure containing numerical solution from ode45 [x, y]

arguments
    a (1,1) double {mustBeNumeric} = 2
    b (1,1) double {mustBeNumeric} = 1
    x0 (1,1) double {mustBeNumeric} = 0
    y0 (1,1) double {mustBeNumeric} = 0
    x_span (1,2) double {mustBeNumeric} = [0, 1.5]
end

%% 1. Analytical Solution Derivation via Separation of Variables
% Mathematical Derivation:
% dy/(a*y + b) = dx
% (1/a) * ln|a*y + b| = x + C1
% a*y + b = C * exp(a*x)  ->  y(x) = C_const * exp(a*x) - b/a
% Enforcing initial state condition y(x0) = y0:
% C_const = (y0 + b/a) * exp(-a*x0)

if abs(a) < eps
    % Special case: Pure integrator (y' = b)
    y_analytic_func = @(x) b * (x - x0) + y0;
    C_const = NaN;
else
    steady_state_offset = -b / a;
    C_const = (y0 - steady_state_offset) * exp(-a * x0);
    y_analytic_func = @(x) C_const * exp(a * x) + steady_state_offset;
end

%% 2. Numerical Integration via ode45 for Verification
ode_fun = @(x, y) a * y + b;
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[x_eval, y_eval] = ode45(ode_fun, x_span, y0, opts);

sol_num.x = x_eval;
sol_num.y = y_eval;

%% 3. Terminal Diagnostics Output
fprintf('====================================================\n');
fprintf('   First-Order ODE Solver: dy/dx = %.4f*y + %.4f   \n', a, b);
fprintf('====================================================\n');
fprintf('Initial Condition     : y(%.4f) = %.4f\n', x0, y0);
if abs(a) >= eps
    fprintf('General Solution Form : y(x) = C * exp(%.4f * x) + (%.4f)\n', a, -b/a);
    fprintf('Particular Solution   : y(x) = %.4f * exp(%.4f * x) + (%.4f)\n', ...
        C_const, a, -b/a);
else
    fprintf('Particular Solution   : y(x) = %.4f * (x - %.4f) + %.4f\n', b, x0, y0);
end

y_ana_eval = y_analytic_func(sol_num.x);
max_error = max(abs(y_ana_eval - sol_num.y));
fprintf('Verification Max Error (Analytic vs ode45): %e\n', max_error);
fprintf('====================================================\n\n');

%% 4. Automatic Demonstration Plot (Active when run via F5 / no output assigned)
if nargout == 0
    figure('Color', [1 1 1], 'Position', [200, 200, 750, 480]);
    x_dense = linspace(x_span(1), x_span(2), 250);

    plot(x_dense, y_analytic_func(x_dense), 'b-', 'LineWidth', 1.8, ...
        'DisplayName', 'Analytic Closed-Form Solution');
    hold on; grid on; box on;
    plot(sol_num.x, sol_num.y, 'ro', 'MarkerSize', 5, 'LineWidth', 1.2, ...
        'DisplayName', 'ode45 Numerical Verification');

    xlabel('Independent Variable x', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('State Variable y(x)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf("Verification: dy/dx = %.2f y + %.2f, y(%.2f) = %.2f", ...
        a, b, x0, y0), 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'northwest', 'FontSize', 10);
    hold off;
end
end