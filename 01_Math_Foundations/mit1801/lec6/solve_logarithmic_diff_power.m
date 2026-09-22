function sol = solve_logarithmic_diff_power(u_sym, v_sym, eval_point)
% SOLVE_LOGARITHMIC_DIFF_POWER Computes symbolic and numerical derivative
% of generalized power function f(x) = u(x)^v(x) using logarithmic differentiation:
% Formula: f'(x) = f(x) * [ v'(x)*ln(u(x)) + v(x)*(u'(x)/u(x)) ]
%
% Syntax:
%   sol = solve_logarithmic_diff_power(u_sym, v_sym, eval_point)
%
% Inputs:
%   u_sym      - Symbolic base function u(x)
%   v_sym      - Symbolic exponent function v(x)
%   eval_point - Operating scalar point x0 (must ensure u(x0) > 0)
%
% Outputs:
%   sol        - Struct containing symbolic expressions, analytic value, and numerical check

arguments
    u_sym      = []
    v_sym      = []
    eval_point (1,1) double {mustBePositive} = 1.5
end

syms x positive
if isempty(u_sym), u_sym = x; end
if isempty(v_sym), v_sym = x; end % Default: f(x) = x^x

% 1. Symbolic logarithmic differentiation setup: ln(f) = v * ln(u)
f_sym = u_sym^v_sym;
du_dx = diff(u_sym, x);
dv_dx = diff(v_sym, x);

d_ln_f = dv_dx * log(u_sym) + v_sym * (du_dx / u_sym);
df_dx_analytic = simplify(f_sym * d_ln_f);

% 2. Numerical evaluation at operating point
f_val = double(subs(f_sym, x, eval_point));
df_analytic_val = double(subs(df_dx_analytic, x, eval_point));
log_derivative_rate = double(subs(d_ln_f, x, eval_point));

% 3. Numerical verification via central finite difference
h = 1e-6;
f_num = matlabFunction(f_sym);
df_numerical_val = (f_num(eval_point + h) - f_num(eval_point - h)) / (2 * h);
abs_discrepancy = abs(df_analytic_val - df_numerical_val);

% Assemble solution structure
sol.x0 = eval_point;
sol.function_expression = char(f_sym);
sol.symbolic_derivative = char(df_dx_analytic);
sol.function_value = f_val;
sol.relative_rate_of_change = log_derivative_rate; % f'(x0) / f(x0)
sol.analytic_derivative = df_analytic_val;
sol.numerical_derivative = df_numerical_val;
sol.discrepancy = abs_discrepancy;

% Console report
fprintf('=== Logarithmic Differentiation of Power Function ===\n');
fprintf('Function Expression f(x)    : %s\n', sol.function_expression);
fprintf('Operating Point (x0)        : %.4f\n', eval_point);
fprintf('Analytical Derivative f''(x) : %s\n', sol.symbolic_derivative);
fprintf('Function Value f(x0)        : %.6f\n', f_val);
fprintf('Relative Growth Rate f''/f   : %.6f\n', log_derivative_rate);
fprintf('Analytical Derivative Value : %.6f\n', df_analytic_val);
fprintf('Numerical Finite Difference : %.6f (Error: %.2e)\n', ...
    df_numerical_val, abs_discrepancy);
end