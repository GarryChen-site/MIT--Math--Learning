function result = solve_tangent_triangle_area(k_param, x0_val)
% SOLVE_TANGENT_TRIANGLE_AREA Computes the tangent line intercepts and
% enclosing triangle area for hyperbola f(x) = k / x.
%
% Syntax:
%   result = solve_tangent_triangle_area(k_param, x0_val)
%
% Inputs:
%   k_param - Proportional constant (f(x) = k/x, default = 1.0)
%   x0_val  - Tangency point x-coordinate (default = 1.0)
%
% Outputs:
%   result  - Struct containing slope, intercepts, area, and symbolic formula

arguments
    k_param (1,1) double {mustBePositive} = 1.0
    x0_val  (1,1) double {mustBePositive} = 1.0
end

syms x x0 k positive
f(x) = k / x;

% 1. Analytic differentiation
df = diff(f, x);
slope_sym = df(x0);

% 2. Intercepts from tangent equation: y - y0 = f'(x0) * (x - x0)
y0 = f(x0);
syms X_int Y_int
eq_x = 0 - y0 == slope_sym * (X_int - x0);
eq_y = Y_int - y0 == slope_sym * (0 - x0);

x_int_sym = simplify(solve(eq_x, X_int));
y_int_sym = simplify(solve(eq_y, Y_int));

% 3. Enclosed area: S = 0.5 * base * height
area_sym = simplify(0.5 * x_int_sym * y_int_sym);

% 4. Numerical substitution
result.k = k_param;
result.x0 = x0_val;
result.slope = double(subs(slope_sym, [k, x0], [k_param, x0_val]));
result.x_intercept = double(subs(x_int_sym, [k, x0], [k_param, x0_val]));
result.y_intercept = double(subs(y_int_sym, [k, x0], [k_param, x0_val]));
result.area = double(subs(area_sym, [k, x0], [k_param, x0_val]));
result.area_symbolic_expression = char(area_sym);

% Console report
fprintf('=== Tangent Triangle Analysis for f(x) = k/x ===\n');
fprintf('Operating Point: x0 = %.4f, k = %.4f\n', x0_val, k_param);
fprintf('Tangent Slope  : %.4f\n', result.slope);
fprintf('X-Intercept    : %.4f (Theoretical: 2 * x0 = %.4f)\n', result.x_intercept, 2 * x0_val);
fprintf('Y-Intercept    : %.4f (Theoretical: 2 * y0 = %.4f)\n', result.y_intercept, 2 * (k_param / x0_val));
fprintf('Enclosed Area  : %.4f (Symbolic Formula: %s)\n', result.area, result.area_symbolic_expression);
end