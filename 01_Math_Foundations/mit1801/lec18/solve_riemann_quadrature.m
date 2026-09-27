function sol = solve_riemann_quadrature(f_handle, a, b, n_steps, f_analytic_integral)
% SOLVE_RIEMANN_QUADRATURE Evaluates definite integral \int_a^b f(x)dx using
% classical Riemann partition methods: Left, Right, Midpoint, and Trapezoidal rules.
%
% Syntax:
%   sol = solve_riemann_quadrature(f_handle, a, b, n_steps, f_analytic_integral)
%
% Inputs:
%   f_handle            - Function handle for integrand f(x)
%   a                   - Lower integration limit (default = 0.0)
%   b                   - Upper integration limit (must be > a, default = 2.0)
%   n_steps             - Number of uniform subintervals (default = 50)
%   f_analytic_integral - Optional analytic solution or integral handle (default = empty)
%
% Outputs:
%   sol                 - Struct containing estimated integral values and errors

arguments
    f_handle            function_handle = @(x) x.^2
    a                   (1,1) double = 0.0
    b                   (1,1) double = 2.0
    n_steps             (1,1) double {mustBeInteger, mustBePositive} = 50
    f_analytic_integral = []
end

if a >= b
    error('Lower limit a must be strictly less than upper limit b.');
end

dx = (b - a) / n_steps;
x_grid = linspace(a, b, n_steps + 1);

% 1. Left-hand Riemann sum: sum(f(x_i)) * dx, i = 1 to n
x_left = x_grid(1:end-1);
left_sum = sum(f_handle(x_left)) * dx;

% 2. Right-hand Riemann sum: sum(f(x_i)) * dx, i = 2 to n+1
x_right = x_grid(2:end);
right_sum = sum(f_handle(x_right)) * dx;

% 3. Midpoint Riemann sum: sum(f(x_mid)) * dx
x_mid = 0.5 * (x_grid(1:end-1) + x_grid(2:end));
midpoint_sum = sum(f_handle(x_mid)) * dx;

% 4. Trapezoidal rule: 0.5 * (Left + Right)
trapezoid_sum = 0.5 * (left_sum + right_sum);

% 5. Determine ground truth reference
if isempty(f_analytic_integral)
    % High-precision adaptive quadrature benchmark
    exact_val = integral(f_handle, a, b, 'RelTol', 1e-12, 'AbsTol', 1e-14);
elseif isa(f_analytic_integral, 'function_handle')
    exact_val = f_analytic_integral(b) - f_analytic_integral(a);
else
    exact_val = double(f_analytic_integral);
end

% Error calculations
err_left = abs(left_sum - exact_val);
err_right = abs(right_sum - exact_val);
err_midpoint = abs(midpoint_sum - exact_val);
err_trapezoid = abs(trapezoid_sum - exact_val);

% Assemble solution structure
sol.interval = [a, b];
sol.n_partitions = n_steps;
sol.step_size_dx = dx;
sol.reference_integral = exact_val;
sol.left_riemann_sum = left_sum;
sol.right_riemann_sum = right_sum;
sol.midpoint_sum = midpoint_sum;
sol.trapezoidal_sum = trapezoid_sum;
sol.errors.left = err_left;
sol.errors.right = err_right;
sol.errors.midpoint = err_midpoint;
sol.errors.trapezoid = err_trapezoid;

% Console report
fprintf('=== Riemann Numerical Quadrature Evaluation ===\n');
fprintf('Integration Interval [a, b]: [%.4f, %.4f] (\\Delta x = %.6f, n = %d)\n', ...
    a, b, dx, n_steps);
fprintf('Reference Integral Value   : %.10f\n', exact_val);
fprintf('Left-Hand Sum (Euler-Fwd)  : %.10f (Error: %.2e)\n', left_sum, err_left);
fprintf('Right-Hand Sum (Euler-Bwd) : %.10f (Error: %.2e)\n', right_sum, err_right);
fprintf('Midpoint Rule (Order 2)    : %.10f (Error: %.2e)\n', midpoint_sum, err_midpoint);
fprintf('Trapezoidal Rule (Tustin)  : %.10f (Error: %.2e)\n', trapezoid_sum, err_trapezoid);
end