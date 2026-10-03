function sol = solve_gaussian_series_integral(upper_limit_x, truncation_terms_N)
% SOLVE_GAUSSIAN_SERIES_INTEGRAL Evaluates the non-elementary Gaussian integral:
% F(x) = \int_0^x e^{-t^2} dt using term-by-term power series integration (Lecture 38).
% Formula: F_N(x) = \sum_{n=0}^N [ (-1)^n * x^(2n+1) ] / [ n! * (2n+1) ]
%
% Syntax:
%   sol = solve_gaussian_series_integral(upper_limit_x, truncation_terms_N)
%
% Inputs:
%   upper_limit_x      - Real evaluation limit coordinate x (default = 1.0)
%   truncation_terms_N - Number of polynomial expansion orders N (default = 5)
%
% Outputs:
%   sol                - Struct containing series value, exact numerical value, and error

arguments
    upper_limit_x      (1,1) double = 1.0
    truncation_terms_N (1,1) double {mustBeInteger, mustBeNonnegative} = 5
end

x_val = upper_limit_x;

% 1. Evaluate Term-by-Term Maclaurin Series Sum
term_values = zeros(truncation_terms_N + 1, 1);
partial_sums = zeros(truncation_terms_N + 1, 1);
running_sum = 0.0;

for n = 0:truncation_terms_N
    % General term: (-1)^n * x^(2n+1) / ( n! * (2n+1) )
    term_val = ((-1)^n * (x_val^(2*n + 1))) / (factorial(n) * (2*n + 1));
    running_sum = running_sum + term_val;

    term_values(n+1) = term_val;
    partial_sums(n+1) = running_sum;
end

series_estimate = running_sum;

% 2. High-Precision Numerical Benchmark via Adaptive Quadrature
f_gaussian = @(t) exp(-t.^2);
exact_integral_val = integral(f_gaussian, 0, x_val, 'RelTol', 1e-12, 'AbsTol', 1e-14);

discrepancy = abs(exact_integral_val - series_estimate);
rel_error_pct = (discrepancy / abs(exact_integral_val)) * 100.0;

% Theoretical Alternating Series Remainder Bound: |R_N| <= |term_{N+1}|
theoretical_tail_bound = abs(((-1)^(truncation_terms_N + 1) * (x_val^(2*truncation_terms_N + 3))) / ...
    (factorial(truncation_terms_N + 1) * (2*truncation_terms_N + 3)));

% Assemble solution structure
sol.upper_limit = x_val;
sol.truncation_orders = truncation_terms_N;
sol.series_terms = term_values;
sol.partial_sums = partial_sums;
sol.series_integral_value = series_estimate;
sol.benchmark_numerical_value = exact_integral_val;
sol.absolute_error = discrepancy;
sol.relative_error_pct = rel_error_pct;
sol.theoretical_remainder_bound = theoretical_tail_bound;

% Console report
fprintf('=== Gaussian Integral via Term-by-Term Series Integration ===\n');
fprintf('Integration Limit (x)       : %.4f (Domain: [0, %.4f])\n', x_val, x_val);
fprintf('Truncation Expansion (N)   : %d Orders (Total %d Odd Polynomial Terms)\n', ...
    truncation_terms_N, truncation_terms_N + 1);
fprintf('Series Approximation Value : %.10f\n', series_estimate);
fprintf('Adaptive Benchmark Quadrature: %.10f\n', exact_integral_val);
fprintf('Absolute Discrepancy Error  : %.2e (Rel: %.4f%%)\n', discrepancy, rel_error_pct);
fprintf('Alternating Series Tail Bound: %.2e (Strictly Bounded: %s)\n', ...
    theoretical_tail_bound, mat2str(discrepancy <= theoretical_tail_bound));
end