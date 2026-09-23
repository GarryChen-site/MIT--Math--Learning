function report = solve_transcendental_peak_kinematics(k_scale, n_power)
% SOLVE_TRANSCENDENTAL_PEAK_KINEMATICS Analyzes transcendental envelope functions
% of the form f(x) = ln(k*x) / x^n for x > 0 (Lecture 10 Example 4 generalization).
% Solves exact peak location, critical value, and horizontal limiting behaviors.
%
% Syntax:
%   report = solve_transcendental_peak_kinematics(k_scale, n_power)
%
% Inputs:
%   k_scale - Positive scaling factor inside logarithm (default = 1.0)
%   n_power - Exponent in denominator (default = 1.0)
%
% Outputs:
%   report  - Struct containing critical coordinates, root, and asymptotic limit

arguments
    k_scale (1,1) double {mustBePositive} = 1.0
    n_power (1,1) double {mustBePositive} = 1.0
end

% 1. Root location: ln(k*x) = 0 => x_zero = 1/k
x_zero = 1.0 / k_scale;

% 2. Derivative derivation:
% f'(x) = [ (1/x)*x^n - n*x^(n-1)*ln(k*x) ] / x^(2n) = [ 1 - n*ln(k*x) ] / x^(n+1)
% f'(x) = 0 => ln(k*x) = 1/n => k*x = exp(1/n) => x_crit = exp(1/n) / k
x_crit = exp(1.0 / n_power) / k_scale;
y_crit = (log(k_scale * x_crit)) / (x_crit^n_power); % Exact value: (1/n) / (e/k^n)

% 3. Inflection point derivation:
% f''(x) = 0 => -n - (n+1)(1 - n*ln(k*x)) = 0
% -(2n+1) + n(n+1)*ln(k*x) = 0 => ln(k*x) = (2n+1)/(n*(n+1))
x_inflect = exp((2.0 * n_power + 1.0) / (n_power * (n_power + 1.0))) / k_scale;
y_inflect = (log(k_scale * x_inflect)) / (x_inflect^n_power);

% Assemble reporting structure
report.k_scale = k_scale;
report.n_power = n_power;
report.root_x = x_zero;
report.critical_point_x = x_crit;
report.maximum_value_y = y_crit;
report.inflection_point_x = x_inflect;
report.inflection_point_y = y_inflect;
report.limit_as_x_to_zero = -Inf;
report.limit_as_x_to_infinity = 0.0;

% Console output
fprintf('=== Transcendental Curve Kinematics: ln(%.2f*x) / x^%.2f ===\n', k_scale, n_power);
fprintf('Real Root (f(x) = 0)       : x = %.4f\n', x_zero);
fprintf('Global Maximum Location    : x = %.4f (Theory: e^(1/%.1f)/%.2f)\n', ...
    x_crit, n_power, k_scale);
fprintf('Global Maximum Amplitude   : y = %.4f\n', y_crit);
fprintf('Inflection Transition Point: (x = %.4f, y = %.4f)\n', x_inflect, y_inflect);
fprintf('Limiting Behaviors         : x -> 0+ => -Inf | x -> +Inf => 0.0000\n');
end