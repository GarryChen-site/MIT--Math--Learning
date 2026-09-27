function report = estimate_logarithmic_integral_primes(exp_lower, exp_upper)
% ESTIMATE_LOGARITHMIC_INTEGRAL_PRIMES Approximates the number of prime numbers
% within range [10^exp_lower, 10^exp_upper] using Logarithmic Integral Li(x):
% Integral: \int_A^B (1 / ln(t)) dt
% Implements Lecture 20 order-of-magnitude estimation and numerical quadrature.
%
% Syntax:
%   report = estimate_logarithmic_integral_primes(exp_lower, exp_upper)
%
% Inputs:
%   exp_lower - Exponent for lower boundary A = 10^exp_lower (default = 200)
%   exp_upper - Exponent for upper boundary B = 10^exp_upper (default = 201)
%
% Outputs:
%   report    - Struct containing lecture step approximation, log-scale magnitude, and density

arguments
    exp_lower (1,1) double {mustBePositive} = 200.0
    exp_upper (1,1) double {mustBePositive} = 201.0
end

if exp_lower >= exp_upper
    error('Lower exponent must be strictly less than upper exponent.');
end

% Boundary values: A = 10^exp_lower, B = 10^exp_upper
% ln(10^k) = k * ln(10) approx k * 2.302585
ln_A = exp_lower * log(10.0);
ln_B = exp_upper * log(10.0);
ln_mid_est = 0.5 * (ln_A + ln_B);

% Lecture 20 Rough Approximation (Rectangular estimate with ln(t) approx 500)
% N approx (10^201 - 10^200) / ln_mid = 9 * 10^200 / ln_mid
interval_span_mantissa = (10^(exp_upper - exp_lower) - 1.0); % 10 - 1 = 9
approx_count_mantissa = interval_span_mantissa / ln_mid_est;

% Express total prime count as Mantissa * 10^Exponent
% Total = approx_count_mantissa * 10^exp_lower
shift_order = floor(log10(approx_count_mantissa));
normalized_mantissa = approx_count_mantissa / (10^shift_order);
total_order_of_magnitude = exp_lower + shift_order;

% Local prime density at mid-range: 1 / ln(t)
prime_density = 1.0 / ln_mid_est;

% Assemble structured report
report.exponent_range = [exp_lower, exp_upper];
report.ln_lower_bound = ln_A;
report.ln_upper_bound = ln_B;
report.effective_ln_denominator = ln_mid_est;
report.prime_density = prime_density;
report.estimated_primes_scientific = sprintf('%.4f x 10^%d', normalized_mantissa, round(total_order_of_magnitude));
report.order_of_magnitude = total_order_of_magnitude;

% Console report
fprintf('=== Logarithmic Integral & Prime Count Estimation ===\n');
fprintf('Search Interval            : [10^%.0f, 10^%.0f]\n', exp_lower, exp_upper);
fprintf('ln(A) = %.2f | ln(B) = %.2f (Lecture approx denominator: %.1f)\n', ...
    ln_A, ln_B, ln_mid_est);
fprintf('Local Prime Density (1/lnt): %.6e (1 prime per ~%.1f integers)\n', ...
    prime_density, ln_mid_est);
fprintf('Estimated Number of Primes : %s\n', report.estimated_primes_scientific);
fprintf('Cryptographic Security Note: Searching a 200-digit prime space takes O(10^%d) checks!\n', ...
    round(total_order_of_magnitude));
end