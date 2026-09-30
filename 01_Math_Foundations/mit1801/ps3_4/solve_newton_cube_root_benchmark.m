function report = solve_newton_cube_root_benchmark(target_val, x0_init, max_iters)
% SOLVE_NEWTON_CUBE_ROOT_BENCHMARK Solves cube root x = target_val^(1/3) using
% Newton-Raphson iteration on f(x) = x^3 - target_val = 0 (PS3 Part II Problem 5).
% Formula: x_{k+1} = (2/3)*x_k + target_val / (3*x_k^2)
% Compares quadratic Taylor perturbation: 9^(1/3) = 2*(1 + 1/8)^(1/3).
%
% Syntax:
%   report = solve_newton_cube_root_benchmark(target_val, x0_init, max_iters)

arguments
    target_val (1,1) double {mustBePositive} = 9.0
    x0_init    (1,1) double {mustBePositive} = 2.0
    max_iters  (1,1) double {mustBeInteger, mustBePositive} = 6
end

exact_root = target_val^(1.0 / 3.0); % Ground truth

% 1. Quadratic Taylor Series Expansion: 9^(1/3) = 8^(1/3) * (1 + 1/8)^(1/3)
% (1 + u)^(1/3) approx 1 + (1/3)*u - (1/9)*u^2
u_pert = (target_val - 8.0) / 8.0; % 1/8 for target 9
taylor_quad = 2.0 * (1.0 + (1.0 / 3.0) * u_pert - (1.0 / 9.0) * (u_pert^2));
taylor_error = abs(taylor_quad - exact_root);

% 2. Newton Iteration Execution
x_hist = zeros(max_iters + 1, 1);
err_hist = zeros(max_iters + 1, 1);

x_curr = x0_init;
x_hist(1) = x_curr;
err_hist(1) = abs(x_curr - exact_root);

for k = 1:max_iters
    x_next = (2.0 / 3.0) * x_curr + target_val / (3.0 * (x_curr^2));
    x_hist(k+1) = x_next;
    err_hist(k+1) = abs(x_next - exact_root);
    x_curr = x_next;
end

% Assemble report structure
report.target_number = target_val;
report.exact_root = exact_root;
report.taylor_quadratic_approx = taylor_quad;
report.taylor_absolute_error = taylor_error;
report.iteration_values = x_hist;
report.iteration_errors = err_hist;

% Console report
fprintf('=== PS3 Problem 5: Newton-Raphson vs. Taylor Cube Root ===\n');
fprintf('Target Base Number          : %.4f (Exact Root = %.12f)\n', target_val, exact_root);
fprintf('2nd-Order Taylor Approx     : %.12f (Error: %.2e)\n', taylor_quad, taylor_error);
fprintf('Newton-Raphson Iterations (Starting x0 = %.1f):\n', x0_init);
for k = 1:length(x_hist)
    fprintf('  Iter %d: x = %16.12f | Residual = %11.2e\n', ...
        k-1, x_hist(k), err_hist(k));
end
end