function sol = solve_differential_approximation(f_sym, base_a, target_x)
% SOLVE_DIFFERENTIAL_APPROXIMATION Estimates f(x) at target_x = base_a + dx
% using differential calculus: f(x) approx f(a) + dy, where dy = f'(a)*dx.
%
% Syntax:
%   sol = solve_differential_approximation(f_sym, base_a, target_x)
%
% Inputs:
%   f_sym    - Symbolic function f(x) (default: x^(1/3))
%   base_a   - Known basepoint with simple evaluation (default: 64.0)
%   target_x - Target evaluation coordinate (default: 65.0)
%
% Outputs:
%   sol      - Struct containing analytic derivative, differential dy, approx vs exact values

arguments
    f_sym    = []
    base_a   (1,1) double = 64.0
    target_x (1,1) double = 65.0
end

syms x
if isempty(f_sym), f_sym = x^(1/3); end

% 1. Symbolic differentiation for differential dy = f'(a) * dx
df_sym = diff(f_sym, x);

dx = target_x - base_a;
f_a = double(subs(f_sym, x, base_a));
df_a = double(subs(df_sym, x, base_a));

% Differential increment
dy = df_a * dx;
approx_val = f_a + dy;

% Exact numeric evaluation
exact_val = double(subs(f_sym, x, target_x));
abs_err = abs(exact_val - approx_val);
rel_err_pct = (abs_err / abs(exact_val)) * 100.0;

% Assemble structured output
sol.function_expression = char(f_sym);
sol.derivative_expression = char(df_sym);
sol.basepoint_a = base_a;
sol.target_x = target_x;
sol.differential_dx = dx;
sol.f_base_a = f_a;
sol.derivative_at_a = df_a;
sol.differential_dy = dy;
sol.approximated_value = approx_val;
sol.exact_value = exact_val;
sol.absolute_error = abs_err;
sol.relative_error_percent = rel_err_pct;

% Console report
fprintf('=== Differential Linear Approximation Analysis ===\n');
fprintf('Target Function f(x)       : %s\n', sol.function_expression);
fprintf('Basepoint (a)              : %.4f => f(a) = %.6f\n', base_a, f_a);
fprintf('Target Evaluation (x)      : %.4f (dx = %+.4f)\n', target_x, dx);
fprintf('Derivative at Base f''(a)   : %.6f\n', df_a);
fprintf('Differential Increment (dy): %.6f\n', dy);
fprintf('Approximated Value f(a)+dy : %.8f\n', approx_val);
fprintf('Exact Function Value f(x)  : %.8f\n', exact_val);
fprintf('Approximation Error        : Absolute = %.2e | Relative = %.4f%%\n', ...
    abs_err, rel_err_pct);
end