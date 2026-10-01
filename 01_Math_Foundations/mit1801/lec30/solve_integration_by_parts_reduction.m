function report = solve_integration_by_parts_reduction(lambda_rate, target_order_n, upper_limit_T)
% SOLVE_INTEGRATION_BY_PARTS_REDUCTION Solves the definite recurrence integral:
% I_n = \int_0^T x^n * exp(lambda * x) dx using integration by parts reduction:
% Formula: I_n = (1 / lambda) * [ T^n * exp(lambda*T) - n * I_{n-1} ], with I_0 = (exp(lambda*T) - 1) / lambda.
% Validates analytical recurrence against adaptive numerical quadrature.
%
% Syntax:
%   report = solve_integration_by_parts_reduction(lambda_rate, target_order_n, upper_limit_T)
%
% Inputs:
%   lambda_rate    - Real exponential decay/growth factor lambda ~= 0 (default = -1.0)
%   target_order_n - Maximum integer power n (default = 4)
%   upper_limit_T  - Upper limit T > 0 (default = 2.0)
%
% Outputs:
%   report         - Struct containing recursion history, exact values, and truncation residuals

arguments
    lambda_rate    (1,1) double {mustBeNonzero} = -1.0
    target_order_n (1,1) double {mustBeInteger, mustBeNonnegative} = 4
    upper_limit_T  (1,1) double {mustBePositive} = 2.0
end

I_recur = zeros(target_order_n + 1, 1);
I_num   = zeros(target_order_n + 1, 1);
err_log = zeros(target_order_n + 1, 1);

% Base case n = 0: \int_0^T exp(lambda*x) dx = (exp(lambda*T) - 1) / lambda
exp_term = exp(lambda_rate * upper_limit_T);
I_recur(1) = (exp_term - 1.0) / lambda_rate;
I_num(1)   = integral(@(x) exp(lambda_rate * x), 0, upper_limit_T);
err_log(1) = abs(I_recur(1) - I_num(1));

% Forward Reduction Loop
for k = 1:target_order_n
    % I_k = (1 / lambda) * [ T^k * e^(lambda*T) - k * I_{k-1} ]
    boundary_uv = (upper_limit_T^k) * exp_term;
    I_recur(k+1) = (boundary_uv - k * I_recur(k)) / lambda_rate;

    % Numerical Benchmark
    I_num(k+1) = integral(@(x) (x.^k) .* exp(lambda_rate * x), 0, upper_limit_T, ...
        'RelTol', 1e-12, 'AbsTol', 1e-14);
    err_log(k+1) = abs(I_recur(k+1) - I_num(k+1));
end

% Assemble structured report
report.lambda_factor = lambda_rate;
report.max_order = target_order_n;
report.upper_limit = upper_limit_T;
report.orders = (0:target_order_n)';
report.recurrence_integral_values = I_recur;
report.numerical_integral_values = I_num;
report.discrepancies = err_log;
report.is_numerically_verified = all(err_log < 1e-10);

% Console report
fprintf('=== Integration by Parts: Recurrence Reduction Solver ===\n');
fprintf('Model Formulation          : I_n = \\int_0^{%.2f} x^n \\cdot e^{%.2fx} dx\n', ...
    upper_limit_T, lambda_rate);
fprintf('Target Polynomial Power    : n = %d\n', target_order_n);
fprintf('Marching Reduction Table:\n');
for k = 0:target_order_n
    fprintf('  Order n = %d: Analytic = %14.10f | Quadrature = %14.10f | Residual = %.2e\n', ...
        k, I_recur(k+1), I_num(k+1), err_log(k+1));
end
fprintf('Algorithmic Stability Check: %s (Max Error: %.2e)\n', ...
    mat2str(report.is_numerically_verified), max(err_log));
end