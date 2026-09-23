function sol = analyze_curve_extrema_and_inflections(f_sym, x_interval)
% ANALYZE_CURVE_EXTREMA_AND_INFLECTIONS Fully analyzes an analytic function:
% calculates real roots (zeros), critical points (f'=0), classifies local extrema
% via the 2nd derivative test (f''), and determines inflection points (f''=0).
%
% Syntax:
%   sol = analyze_curve_extrema_and_inflections(f_sym, x_interval)
%
% Inputs:
%   f_sym      - Symbolic scalar expression in variable 'x' (default: 3*x - x^3)
%   x_interval - Real search interval [x_min, x_max] (default: [-3, 3])
%
% Outputs:
%   sol        - Struct containing real roots, critical point classification, and inflections

arguments
    f_sym      = []
    x_interval (1,2) double = [-3.0, 3.0]
end

syms x real
if isempty(f_sym), f_sym = 3*x - x^3; end

% 1. Successive symbolic derivatives
df_sym  = diff(f_sym, x);
d2f_sym = diff(df_sym, x);

% 2. Find real roots (zeros of f(x) = 0)
roots_all = double(solve(f_sym == 0, x, 'Real', true));
valid_roots = roots_all(roots_all >= x_interval(1) & roots_all <= x_interval(2));

% 3. Find critical points (f'(x) = 0)
crit_all = double(solve(df_sym == 0, x, 'Real', true));
crit_in_range = crit_all(crit_all >= x_interval(1) & crit_all <= x_interval(2));

num_crit = length(crit_in_range);
crit_classification = cell(num_crit, 1);
crit_y_vals = zeros(num_crit, 1);
d2f_at_crit = zeros(num_crit, 1);

for i = 1:num_crit
    xc = crit_in_range(i);
    yc = double(subs(f_sym, x, xc));
    d2fc = double(subs(d2f_sym, x, xc));

    crit_y_vals(i) = yc;
    d2f_at_crit(i) = d2fc;

    % Second Derivative Test Classification
    if d2fc > 1e-9
        crit_classification{i} = 'Local Minimum (Concave Up)';
    elseif d2fc < -1e-9
        crit_classification{i} = 'Local Maximum (Concave Down)';
    else
        crit_classification{i} = 'Inconclusive / Neutral Point (f''''=0)';
    end
end

% 4. Find inflection points (f''(x) = 0)
inflect_all = double(solve(d2f_sym == 0, x, 'Real', true));
inflect_in_range = inflect_all(inflect_all >= x_interval(1) & inflect_all <= x_interval(2));
inflect_y_vals = zeros(length(inflect_in_range), 1);
for j = 1:length(inflect_in_range)
    inflect_y_vals(j) = double(subs(f_sym, x, inflect_in_range(j)));
end

% Assemble structured output
sol.function_expression = char(f_sym);
sol.first_derivative = char(simplify(df_sym));
sol.second_derivative = char(simplify(d2f_sym));
sol.search_interval = x_interval;
sol.real_zeros = sort(valid_roots);
sol.critical_points_x = crit_in_range;
sol.critical_points_y = crit_y_vals;
sol.second_derivative_values = d2f_at_crit;
sol.extrema_types = crit_classification;
sol.inflection_points_x = inflect_in_range;
sol.inflection_points_y = inflect_y_vals;

% Console report
fprintf('=== Curve Analysis & Extrema Classification ===\n');
fprintf('Function f(x)             : %s\n', sol.function_expression);
fprintf('f''(x)                     : %s\n', sol.first_derivative);
fprintf('f''''(x)                    : %s\n', sol.second_derivative);
fprintf('Search Domain             : [%.2f, %.2f]\n', x_interval(1), x_interval(2));
fprintf('Real Zeros (f(x) = 0)     : %s\n', mat2str(sol.real_zeros', 4));
fprintf('Identified Critical Points:\n');
for i = 1:num_crit
    fprintf('  Point %d: (x = %7.4f, y = %7.4f) | f'''' = %7.4f -> %s\n', ...
        i, crit_in_range(i), crit_y_vals(i), d2f_at_crit(i), crit_classification{i});
end
fprintf('Inflection Points (f''''=0) :\n');
for j = 1:length(inflect_in_range)
    fprintf('  Inflection %d: (x = %7.4f, y = %7.4f)\n', ...
        j, inflect_in_range(j), inflect_y_vals(j));
end
end