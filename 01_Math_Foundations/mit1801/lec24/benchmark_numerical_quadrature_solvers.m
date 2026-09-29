function benchmark_table = benchmark_numerical_quadrature_solvers(f_handle, a, b, exact_val)
% BENCHMARK_NUMERICAL_QUADRATURE_SOLVERS Benchmarks Left Riemann, Trapezoidal,
% and Simpson's methods across powers-of-two refinements and extracts empirical
% numerical orders of convergence: Error approx C * (Delta x)^p.
%
% Syntax:
%   benchmark_table = benchmark_numerical_quadrature_solvers(f_handle, a, b, exact_val)
%
% Inputs:
%   f_handle  - Function handle for integrand f(x)
%   a         - Lower limit (default = 0.0)
%   b         - Upper limit (default = 1.0)
%   exact_val - Ground truth value (default = pi/4)
%
% Outputs:
%   benchmark_table - Table containing steps, values, errors, and convergence orders

arguments
    f_handle  function_handle = @(x) 1.0 ./ (1.0 + x.^2)
    a         (1,1) double = 0.0
    b         (1,1) double = 1.0
    exact_val (1,1) double = pi / 4.0
end

n_cases = [2, 4, 8, 16, 32, 64]';
n_rows = length(n_cases);

dx_col    = zeros(n_rows, 1);
err_riem  = zeros(n_rows, 1);
err_trap  = zeros(n_rows, 1);
err_simp  = zeros(n_rows, 1);

for k = 1:n_rows
    n = n_cases(k);
    h = (b - a) / n;
    dx_col(k) = h;

    x_mesh = linspace(a, b, n + 1);
    y_mesh = f_handle(x_mesh);

    % Left Riemann
    val_r = h * sum(y_mesh(1:end-1));
    err_riem(k) = abs(val_r - exact_val);

    % Trapezoidal
    val_t = h * (0.5 * y_mesh(1) + sum(y_mesh(2:end-1)) + 0.5 * y_mesh(end));
    err_trap(k) = abs(val_t - exact_val);

    % Simpson's
    val_s = (h / 3.0) * (y_mesh(1) + 4 * sum(y_mesh(2:2:end-1)) + 2 * sum(y_mesh(3:2:end-2)) + y_mesh(end));
    err_simp(k) = abs(val_s - exact_val);
end

% Empirical convergence order: p = log(Err_{k-1} / Err_k) / log(dx_{k-1} / dx_k)
order_riem = [NaN; diff(log(err_riem)) ./ diff(log(dx_col))];
order_trap = [NaN; diff(log(err_trap)) ./ diff(log(dx_col))];
order_simp = [NaN; diff(log(err_simp)) ./ diff(log(dx_col))];

benchmark_table = table(n_cases, dx_col, err_riem, order_riem, ...
    err_trap, order_trap, err_simp, order_simp, ...
    'VariableNames', {'Subdivisions_n', 'Delta_x', ...
    'Error_Riemann', 'Order_Riemann', ...
    'Error_Trapezoid', 'Order_Trapezoid', ...
    'Error_Simpson', 'Order_Simpson'});

% Console report
fprintf('=== Multi-Solver Numerical Quadrature Benchmark ===\n');
fprintf('Benchmark Integrand: f(x) = 1/(1+x^2) on [%.2f, %.2f] (Exact: %.6f)\n\n', a, b, exact_val);
disp(benchmark_table);
fprintf('Theoretical Expected Orders : Riemann = O(dx^1) | Trapezoidal = O(dx^2) | Simpson = O(dx^4)\n');
end