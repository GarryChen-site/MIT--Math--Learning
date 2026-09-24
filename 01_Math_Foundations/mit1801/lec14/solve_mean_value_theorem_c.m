function sol = solve_mean_value_theorem_c(f_sym, interval_ab)
% SOLVE_MEAN_VALUE_THEOREM_C Verifies the Mean Value Theorem (MVT) for function f(x)
% over interval [a, b] and computes all intermediate coordinates c in (a, b)
% such that: f'(c) = [f(b) - f(a)] / (b - a).
%
% Syntax:
%   sol = solve_mean_value_theorem_c(f_sym, interval_ab)
%
% Inputs:
%   f_sym       - Symbolic expression in variable 'x' (default: x^3 - 3*x^2 + 2)
%   interval_ab - Two-element vector [a, b] with a < b (default: [0.5, 3.2])
%
% Outputs:
%   sol         - Struct containing secant slope, evaluated c values, and tangent residuals

arguments
    f_sym       = []
    interval_ab (1,2) double = [0.5, 3.2]
end

syms x real
if isempty(f_sym), f_sym = x^3 - 3*x^2 + 2; end

a = interval_ab(1);
b = interval_ab(2);

if a >= b
    error('Invalid interval: Lower bound a must be strictly less than upper bound b.');
end

% 1. Evaluate boundary values and secant slope
fa = double(subs(f_sym, x, a));
fb = double(subs(f_sym, x, b));
secant_slope = (fb - fa) / (b - a);

% 2. Symbolic derivative and solve for critical points c
df_sym = diff(f_sym, x);
mvt_equation = (df_sym == secant_slope);

% Solve equation analytically or numerically
c_candidates = double(solve(mvt_equation, x, 'Real', true));

% 3. Filter points strictly inside the open interval (a, b)
valid_c_points = c_candidates(c_candidates > a & c_candidates < b);
num_c = length(valid_c_points);

if num_c == 0
    warning('No valid interior point c found. Verify function differentiability on (a, b).');
end

% Verify tangent slope residuals at found points
residuals = zeros(num_c, 1);
for k = 1:num_c
    df_at_c = double(subs(df_sym, x, valid_c_points(k)));
    residuals(k) = abs(df_at_c - secant_slope);
end

% Assemble structured output
sol.function_expression = char(f_sym);
sol.derivative_expression = char(df_sym);
sol.interval = [a, b];
sol.f_a = fa;
sol.f_b = fb;
sol.secant_slope = secant_slope;
sol.intermediate_points_c = valid_c_points;
sol.num_valid_points = num_c;
sol.tangent_discrepancies = residuals;

% Console report
fprintf('=== Mean Value Theorem (MVT) Verification ===\n');
fprintf('Function f(x)             : %s\n', sol.function_expression);
fprintf('Derivative f''(x)          : %s\n', sol.derivative_expression);
fprintf('Interval [a, b]           : [%.4f, %.4f]\n', a, b);
fprintf('Endpoints Evaluation      : f(a) = %.4f | f(b) = %.4f\n', fa, fb);
fprintf('Secant Line Slope         : %.6f\n', secant_slope);
fprintf('Valid Interior Points c in (a, b):\n');
for k = 1:num_c
    fprintf('  Point %d: c = %.6f | f''(c) = %.6f | Error = %.2e\n', ...
        k, valid_c_points(k), secant_slope + residuals(k), residuals(k));
end
end