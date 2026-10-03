function report = evaluate_series_convergence_integral_test(p_exponent, N_partial_terms)
% EVALUATE_SERIES_CONVERGENCE_INTEGRAL_TEST Evaluates convergence of p-series:
% S = \sum_{n=1}^\infty (1 / n^p) using the Integral Comparison Test (Lecture 36).
%
% Formulas:
%   Integral Equivalent : \int_1^\infty x^{-p} dx (Converges for p > 1, Diverges for p <= 1)
%   Integral Bounds     : \int_1^{N+1} x^{-p} dx <= S_N <= 1 + \int_1^N x^{-p} dx
%   Known Benchmark     : S = pi^2 / 6 approx 1.644934 (for p = 2)
%
% Syntax:
%   report = evaluate_series_convergence_integral_test(p_exponent, N_partial_terms)
%
% Inputs:
%   p_exponent      - Exponent power parameter p (default = 2.0 for Basel series)
%   N_partial_terms - Number of finite partial sum terms N (default = 1000)
%
% Outputs:
%   report          - Struct containing convergence verdict, partial sum, and integral bounds

arguments
    p_exponent      (1,1) double {mustBePositive} = 2.0
    N_partial_terms (1,1) double {mustBeInteger, mustBePositive} = 1000
end

p = p_exponent;
N = N_partial_terms;

% 1. Convergence Classification via Integral Test: \int_1^\infty x^(-p) dx
is_convergent = (p > 1.0);

% 2. Finite Partial Sum S_N = \sum_{n=1}^N (1 / n^p)
n_seq = 1:N;
terms_vec = 1.0 ./ (n_seq .^ p);
S_N = sum(terms_vec);

% 3. Integral Comparison Bounds Calculation
if abs(p - 1.0) < 1e-12
    % Harmonic Series: p = 1 => Integral = ln(x)
    lower_bound = log(N + 1.0);
    upper_bound = 1.0 + log(N);
    analytic_infinite_sum = Inf;
else
    % p ~= 1: \int x^(-p) dx = x^(1-p) / (1-p)
    anti_func = @(x) (x .^ (1.0 - p)) / (1.0 - p);

    % Lower integral bound: \int_1^{N+1} f(x)dx
    lower_bound = anti_func(N + 1.0) - anti_func(1.0);

    % Upper integral bound: f(1) + \int_1^N f(x)dx
    upper_bound = 1.0 + (anti_func(N) - anti_func(1.0));

    if is_convergent
        % Analytical limit for p = 2 (Basel Problem: pi^2 / 6)
        if abs(p - 2.0) < 1e-12
            analytic_infinite_sum = (pi^2) / 6.0;
        else
            % High-precision symbolic/numeric approximation
            analytic_infinite_sum = double(zeta(p));
        end
    else
        analytic_infinite_sum = Inf;
    end
end

% Truncation error if convergent
if is_convergent
    truncation_gap = abs(analytic_infinite_sum - S_N);
    % Theoretical tail bound via integral: \int_N^\infty x^(-p) dx = N^(1-p) / (p - 1)
    tail_bound_integral = (N^(1.0 - p)) / (p - 1.0);
else
    truncation_gap = Inf;
    tail_bound_integral = Inf;
end

% Assemble structured report
report.power_p = p;
report.partial_terms_count = N;
report.is_series_convergent = is_convergent;
report.partial_sum_S_N = S_N;
report.lower_integral_bound = lower_bound;
report.upper_integral_bound = upper_bound;
report.is_bounded_verified = (S_N >= lower_bound - 1e-9) && (S_N <= upper_bound + 1e-9);
report.analytic_infinite_sum = analytic_infinite_sum;
report.actual_truncation_error = truncation_gap;
report.theoretical_tail_bound = tail_bound_integral;

% Console report
fprintf('=== Series Convergence & Integral Comparison Test ===\n');
fprintf('Target Series Formulation  : \\sum_{n=1}^\\infty (1 / n^{%.2f})\n', p);
fprintf('Integral Test Verdict      : %s (Criterion: p = %.2f > 1.0)\n', ...
    mat2str(is_convergent), p);
fprintf('Computed Partial Sum S_%d : %.10f\n', N, S_N);
fprintf('Integral Bracket Window    : [%.8f, %.8f] -> Verified: %s\n', ...
    lower_bound, upper_bound, mat2str(report.is_bounded_verified));
if is_convergent
    fprintf('Analytic Sum Limit S_inf   : %.10f\n', analytic_infinite_sum);
    fprintf('Actual Truncation Error    : %.2e (N = %d terms)\n', truncation_gap, N);
    fprintf('Integral Tail Upper Bound  : %.2e (Formula: N^{1-p}/(p-1))\n', tail_bound_integral);
else
    fprintf('Infinite Sum Behavior      : Diverges to +Infinity (Harmonic Rate)\n');
end
end