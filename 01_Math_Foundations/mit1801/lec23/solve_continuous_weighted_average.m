function sol = solve_continuous_weighted_average(f_sym, w_sym, var_bounds)
% SOLVE_CONTINUOUS_WEIGHTED_AVERAGE Computes both uniform and weighted continuous
% averages of function f(x) over interval [a, b] with weight distribution w(x).
% Formulas:
%   Uniform Mean  : \bar{f}_{unif} = \frac{1}{b - a} \int_a^b f(x) dx
%   Weighted Mean : \bar{f}_{weight} = \frac{\int_a^b f(x) w(x) dx}{\int_a^b w(x) dx}
%
% Syntax:
%   sol = solve_continuous_weighted_average(f_sym, w_sym, var_bounds)
%
% Inputs:
%   f_sym      - Symbolic profile function in variable 'x' (default: 100 - 30*x for Cauldron Temp)
%   w_sym      - Symbolic weight function in variable 'x' (default: pi*x for Cauldron Disks)
%   var_bounds - Two-element bounds vector [a, b] (default: [0.0, 1.0])
%
% Outputs:
%   sol        - Struct containing uniform average, weighted average, and integrals

arguments
    f_sym      = []
    w_sym      = []
    var_bounds (1,2) double = [0.0, 1.0]
end

syms x real
if isempty(f_sym), f_sym = 100 - 30*x; end
if isempty(w_sym), w_sym = pi * x; end

a = var_bounds(1);
b = var_bounds(2);

if a >= b
    error('Lower bound a must be strictly less than upper bound b.');
end

span = b - a;

% 1. Uniform Continuous Average: (1 / span) * \int_a^b f(x) dx
int_f_sym = int(f_sym, x, a, b);
unif_mean = double(int_f_sym) / span;

% 2. Weighted Continuous Average: \int f*w dx / \int w dx
int_weighted_num_sym = int(f_sym * w_sym, x, a, b);
int_weight_den_sym = int(w_sym, x, a, b);

weighted_numerator = double(int_weighted_num_sym);
total_weight = double(int_weight_den_sym);

if abs(total_weight) < 1e-14
    error('Total integrated weight is zero; weighted average is undefined.');
end

weighted_mean = weighted_numerator / total_weight;

% Assemble solution structure
sol.variable_bounds = [a, b];
sol.function_expression = char(f_sym);
sol.weight_expression = char(w_sym);
sol.uniform_average = unif_mean;
sol.weighted_average = weighted_mean;
sol.total_weight_integral = total_weight;
sol.weight_bias_shift = weighted_mean - unif_mean;

% Console report
fprintf('=== Continuous & Weighted Average Analysis ===\n');
fprintf('Target Profile f(x)       : %s\n', sol.function_expression);
fprintf('Weighting Metric w(x)     : %s\n', sol.weight_expression);
fprintf('Integration Interval [a, b]: [%.4f, %.4f]\n', a, b);
fprintf('Uniform Continuous Average : %.6f\n', unif_mean);
fprintf('Weighted Space Average     : %.6f\n', weighted_mean);
fprintf('Weight-Induced Bias Shift  : %+.6f\n', sol.weight_bias_shift);
end