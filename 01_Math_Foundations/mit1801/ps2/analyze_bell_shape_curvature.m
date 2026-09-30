function report = analyze_bell_shape_curvature()
% ANALYZE_BELL_SHAPE_CURVATURE Solves critical points, maximum values, and
% inflection points for the benchmark Runge bell-shaped curve: f(x) = 1 / (1 + x^2).
%
% Formulas:
%   f'(x)   = -2*x / (1 + x^2)^2
%   f''(x)  = (6*x^2 - 2) / (1 + x^2)^3
%
% Syntax:
%   report = analyze_bell_shape_curvature()
%
% Outputs:
%   report  - Struct containing critical coordinates, inflection points, and slope extrema

syms x real
f_sym = 1 / (1 + x^2);
df_sym = diff(f_sym, x);
d2f_sym = diff(df_sym, x);

% 1. Critical point: f'(x) = 0 => -2*x = 0 => x = 0
x_crit = double(solve(df_sym == 0, x));
y_crit = double(subs(f_sym, x, x_crit));

% 2. Inflection points: f''(x) = 0 => 6*x^2 - 2 = 0 => x = \pm 1 / sqrt(3)
x_inflect = double(solve(d2f_sym == 0, x));
y_inflect = double(subs(f_sym, x, x_inflect));
slope_at_inflect = double(subs(df_sym, x, x_inflect));

report.function_expression = char(f_sym);
report.first_derivative = char(simplify(df_sym));
report.second_derivative = char(simplify(d2f_sym));
report.critical_point = [x_crit, y_crit];
report.inflection_points_x = x_inflect;
report.inflection_points_y = y_inflect;
report.steepest_slopes = slope_at_inflect;

% Console report
fprintf('=== PS2B Problem 2: Bell-Shaped Curve Curvature Analysis ===\n');
fprintf('Base Function f(x)   : %s\n', report.function_expression);
fprintf('f''(x)                : %s\n', report.first_derivative);
fprintf('f''''(x)               : %s\n', report.second_derivative);
fprintf('Global Maximum       : (x = %.4f, y = %.4f)\n', x_crit, y_crit);
for k = 1:length(x_inflect)
    fprintf('Inflection Point %d   : (x = %+.4f, y = %.4f) | Extreme Slope f''(x) = %+.4f\n', ...
        k, x_inflect(k), y_inflect(k), slope_at_inflect(k));
end
end