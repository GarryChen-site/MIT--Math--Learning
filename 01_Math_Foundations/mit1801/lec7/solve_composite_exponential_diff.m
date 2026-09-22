function result = solve_composite_exponential_diff(u_sym, v_sym, x_eval)
% SOLVE_COMPOSITE_EXPONENTIAL_DIFF Solves composite derivative for functions of
% the form f(x) = exp(u(x) * v(x)) via product and chain rules.
%
% Syntax:
%   result = solve_composite_exponential_diff(u_sym, v_sym, x_eval)
%
% Inputs:
%   u_sym  - Symbolic sub-function u(x) (default = x)
%   v_sym  - Symbolic sub-function v(x) (default = atan(x))
%   x_eval - Point of numerical evaluation (default = 1.0)
%
% Outputs:
%   result - Struct containing symbolic derivative, decomposed terms, and error

arguments
    u_sym  = []
    v_sym  = []
    x_eval (1,1) double = 1.0
end

syms x
if isempty(u_sym), u_sym = x; end
if isempty(v_sym), v_sym = atan(x); end

% Internal product w(x) = u(x) * v(x)
w_sym = u_sym * v_sym;
dw_dx = diff(w_sym, x); % u'*v + u*v'

% Total function f(x) = exp(w(x))
f_sym = exp(w_sym);
df_dx_analytic = simplify(f_sym * dw_dx);

% Numerical substitution
f_val = double(subs(f_sym, x, x_eval));
df_val_analytic = double(subs(df_dx_analytic, x, x_eval));

% High-precision central finite difference
h = 1e-6;
f_handle = matlabFunction(f_sym);
df_val_num = (f_handle(x_eval + h) - f_handle(x_eval - h)) / (2 * h);
abs_err = abs(df_val_analytic - df_val_num);

result.x_eval = x_eval;
result.f_expression = char(f_sym);
result.symbolic_derivative = char(df_dx_analytic);
result.f_value = f_val;
result.analytic_derivative = df_val_analytic;
result.numerical_derivative = df_val_num;
result.validation_discrepancy = abs_err;

% Console output
fprintf('=== Lecture 7 Grand Challenge Solver: d/dx [exp(u*v)] ===\n');
fprintf('Function Expression f(x)   : %s\n', result.f_expression);
fprintf('Symbolic Derivative f''(x)  : %s\n', result.symbolic_derivative);
fprintf('Evaluated at x = %.4f      : f(x) = %.6f\n', x_eval, f_val);
fprintf('Analytic Derivative Value  : %.8f\n', df_val_analytic);
fprintf('Numerical Derivative Value : %.8f (Error: %.2e)\n', df_val_num, abs_err);
end