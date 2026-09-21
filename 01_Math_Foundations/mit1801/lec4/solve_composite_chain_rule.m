function sol = solve_composite_chain_rule(f_sym, g_sym, eval_point)
% SOLVE_COMPOSITE_CHAIN_RULE Computes analytic derivative of composite function
% y = f(g(x)) via Chain Rule and decomposes sub-block sensitivities.
%
% Syntax:
%   sol = solve_composite_chain_rule(f_sym, g_sym, eval_point)
%
% Inputs:
%   f_sym      - Outer function f(u) expressed symbolically with variable 'u'
%   g_sym      - Inner function g(x) expressed symbolically with variable 'x'
%   eval_point - Operating point x0 for numeric evaluation
%
% Outputs:
%   sol        - Struct containing chain-rule terms, analytic value, and validation error

arguments
    f_sym      = []
    g_sym      = []
    eval_point (1,1) double = 1.0
end

syms x u
if isempty(f_sym), f_sym = sin(u); end
if isempty(g_sym), g_sym = x^2; end

% 1. Intermediate state evaluation
u_sym = g_sym;
y_composite = subs(f_sym, u, u_sym);

% 2. Chain rule components
df_du = diff(f_sym, u);
dg_dx = diff(g_sym, x);

% Substitute u = g(x) into df/du
df_du_in_x = subs(df_du, u, u_sym);
dy_dx_chain = df_du_in_x * dg_dx;

% 3. Numerical evaluations at evaluation point
u_val = double(subs(u_sym, x, eval_point));
df_du_val = double(subs(df_du, u, u_val));
dg_dx_val = double(subs(dg_dx, x, eval_point));
dy_dx_val = double(subs(dy_dx_chain, x, eval_point));

% 4. Finite-difference verification (central difference)
h = 1e-6;
f_num = @(val) double(subs(y_composite, x, val));
dy_dx_num = (f_num(eval_point + h) - f_num(eval_point - h)) / (2 * h);
abs_error = abs(dy_dx_val - dy_dx_num);

% Assemble solution structure
sol.x0 = eval_point;
sol.intermediate_u0 = u_val;
sol.outer_derivative_df_du = df_du_val;
sol.inner_derivative_dg_dx = dg_dx_val;
sol.composite_derivative_analytic = dy_dx_val;
sol.composite_derivative_numerical = dy_dx_num;
sol.absolute_discrepancy = abs_error;
sol.symbolic_derivative = char(simplify(dy_dx_chain));

% Console output
fprintf('=== Chain Rule Symbolic & Numerical Decomposition ===\n');
fprintf('Composite Function y(x)    : %s\n', char(y_composite));
fprintf('Outer Block dy/du          : %s\n', char(df_du));
fprintf('Inner Block du/dx          : %s\n', char(dg_dx));
fprintf('Chain Rule Formula dy/dx   : %s\n', sol.symbolic_derivative);
fprintf('Evaluated at x0 = %.4f     : u0 = %.4f\n', eval_point, u_val);
fprintf('Subsystem Sensitivities    : dy/du = %.6f | du/dx = %.6f\n', df_du_val, dg_dx_val);
fprintf('Total Derivative           : %.6f (Numerical: %.6f, Error: %.2e)\n', ...
    dy_dx_val, dy_dx_num, abs_error);
end