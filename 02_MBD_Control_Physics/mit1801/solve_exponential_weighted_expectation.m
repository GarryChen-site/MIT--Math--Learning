function report = solve_exponential_weighted_expectation(half_life_H, payoff_function_handle, finite_horizon_M)
% SOLVE_EXPONENTIAL_WEIGHTED_EXPECTATION Solves the continuous exponentially-weighted
% expected value and functional safety lifetime (MTBF / annuity problem) from Lecture 38:
% Formula: E[f(t)] = [ \int_0^\infty f(t) * e^{-kt} dt ] / [ \int_0^\infty e^{-kt} dt ]
% With rate parameter k = ln(2) / H.
% For linear payoff f(t) = t: E[t] = 1/k = H / ln(2).
%
% Syntax:
%   report = solve_exponential_weighted_expectation(half_life_H, payoff_function_handle, finite_horizon_M)
%
% Inputs:
%   half_life_H            - Half-life H > 0 in seconds or years (default = 2.0)
%   payoff_function_handle - Function handle for payoff/state f(t) (default: @(t) t)
%   finite_horizon_M       - Finite horizon truncation time M (default = 15.0)
%
% Outputs:
%   report                 - Struct containing rate k, MTBF, analytical expectation, and truncated value

arguments
    half_life_H            (1,1) double {mustBePositive} = 2.0
    payoff_function_handle function_handle = @(t) t
    finite_horizon_M       (1,1) double {mustBePositive} = 15.0
end

% 1. Parameter Calibration from Lecture 38: k = ln(2) / H
k_rate = log(2.0) / half_life_H;
tau_mean = 1.0 / k_rate; % Exact MTBF / Expected Lifetime for f(t) = t

% 2. Infinite Horizon Integration Evaluation
% Denominator: \int_0^\infty e^{-kt} dt = 1/k
weight_integral_inf = 1.0 / k_rate;

% Numerator: \int_0^\infty f(t) * e^{-kt} dt
kernel_inf = @(t) payoff_function_handle(t) .* exp(-k_rate * t);
num_integral_inf = integral(kernel_inf, 0, Inf, 'RelTol', 1e-12, 'AbsTol', 1e-14);

expected_value_inf = num_integral_inf / weight_integral_inf;

% 3. Finite Horizon Truncated Integration: [0, M]
weight_integral_M = integral(@(t) exp(-k_rate * t), 0, finite_horizon_M);
num_integral_M = integral(kernel_inf, 0, finite_horizon_M);
expected_value_finite_M = num_integral_M / weight_integral_M;

truncation_residual = abs(expected_value_inf - expected_value_finite_M);

% Assemble structured report
report.half_life_H = half_life_H;
report.decay_constant_k = k_rate;
report.mean_lifetime_mtbf = tau_mean;
report.infinite_expected_value = expected_value_inf;
report.finite_horizon_M = finite_horizon_M;
report.finite_horizon_expected_value = expected_value_finite_M;
report.truncation_residual_error = truncation_residual;

% Console report
fprintf('=== Exponential Weighted Expectation & Lifecycle Solver ===\n');
fprintf('Input Half-Life (H)         : %.4f units\n', half_life_H);
fprintf('Calibrated Decay Factor (k) : %.6f units^-1 (Formula: ln2 / H)\n', k_rate);
fprintf('Theoretical Mean Life (MTBF): %.6f units (Formula: 1 / k = H / ln2)\n', tau_mean);
fprintf('Infinite Horizon Expectation: %.10f\n', expected_value_inf);
fprintf('Finite Truncation at M=%.1f : %.10f\n', finite_horizon_M, expected_value_finite_M);
fprintf('Truncation Discrepancy Error: %.2e\n', truncation_residual);
end