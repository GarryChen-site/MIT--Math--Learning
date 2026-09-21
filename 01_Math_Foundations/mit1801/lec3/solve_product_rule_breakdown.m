function result = solve_product_rule_breakdown(u_func, v_func, x0, dx)
% SOLVE_PRODUCT_RULE_BREAKDOWN Evaluates product rule components and decomposes
% total area change into first-order differentials and second-order cross-term.
%
% Syntax:
%   result = solve_product_rule_breakdown(u_func, v_func, x0, dx)
%
% Inputs:
%   u_func - Function handle for u(x)
%   v_func - Function handle for v(x)
%   x0     - Operating expansion point (default = 1.0)
%   dx     - Perturbation step Delta x (default = 1e-3)
%
% Outputs:
%   result - Struct containing first-order terms, second-order terms, and relative error

arguments
    u_func (1,1) function_handle = @(x) x.^2
    v_func (1,1) function_handle = @(x) sin(x)
    x0     (1,1) double = 1.0
    dx     (1,1) double {mustBePositive} = 1e-3
end

% Base state evaluation
u0 = u_func(x0);
v0 = v_func(x0);
u_pert = u_func(x0 + dx);
v_pert = v_func(x0 + dx);

delta_u = u_pert - u0;
delta_v = v_pert - v0;

% Area terms decomposition
exact_delta_prod = (u_pert * v_pert) - (u0 * v0);
first_order_linear = (u0 * delta_v) + (v0 * delta_u);
second_order_error = delta_u * delta_v;

% Relative contribution of the higher-order infinitesimal
higher_order_ratio = abs(second_order_error / exact_delta_prod);

% Struct assembly
result.x0 = x0;
result.dx = dx;
result.u0 = u0;
result.v0 = v0;
result.delta_u = delta_u;
result.delta_v = delta_v;
result.exact_delta_product = exact_delta_prod;
result.first_order_differential = first_order_linear;
result.second_order_residual = second_order_error;
result.higher_order_ratio = higher_order_ratio;

% Console output
fprintf('=== Product Rule Differential Breakdown ===\n');
fprintf('Expansion Point x0         : %.4f, Step dx: %.2e\n', x0, dx);
fprintf('Exact Total Delta(uv)      : %.8e\n', exact_delta_prod);
fprintf('1st-Order Terms (u*dv+v*du): %.8e\n', first_order_linear);
fprintf('2nd-Order Term  (du*dv)    : %.8e\n', second_order_error);
fprintf('Truncation Error Ratio     : %.4f%%\n', higher_order_ratio * 100);
end