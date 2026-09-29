function sol = solve_simpson_quadrature(f_handle, a, b, num_intervals, exact_reference)
% SOLVE_SIMPSON_QUADRATURE Evaluates definite integral \int_a^b f(x) dx using
% Simpson's 1/3 parabolic rule with automated composite weight synthesis.
% Formula: S_n = (Delta x / 3) * [ y0 + 4*y1 + 2*y2 + 4*y3 + ... + 4*y_{n-1} + yn ]
%
% Syntax:
%   sol = solve_simpson_quadrature(f_handle, a, b, num_intervals, exact_reference)
%
% Inputs:
%   f_handle        - Function handle for integrand f(x)
%   a               - Lower integration boundary (default = 0.0)
%   b               - Upper integration boundary (must be > a, default = 1.0)
%   num_intervals   - Number of intervals n (MUST BE EVEN, default = 2 from lecture)
%   exact_reference - Optional theoretical analytic value (default = pi/4)
%
% Outputs:
%   sol             - Struct containing estimated value, step size, node weights, and error

arguments
    f_handle        function_handle = @(x) 1.0 ./ (1.0 + x.^2)
    a               (1,1) double = 0.0
    b               (1,1) double = 1.0
    num_intervals   (1,1) double {mustBeInteger, mustBePositive} = 2
    exact_reference = []
end

if a >= b
    error('Lower limit a must be strictly less than upper limit b.');
end

% Constraint check: n must be an even integer for Simpson's 1/3 rule
if mod(num_intervals, 2) ~= 0
    error('Simpson''s Rule requires an EVEN number of intervals n (received n = %d).', num_intervals);
end

dx = (b - a) / num_intervals;
x_grid = linspace(a, b, num_intervals + 1);
y_vals = f_handle(x_grid);

% Construct Simpson weight multipliers vector: [1, 4, 2, 4, 2, ..., 4, 1]
weights = 2 * ones(1, num_intervals + 1);
weights(1) = 1;
weights(end) = 1;
weights(2:2:end-1) = 4; % Odd indices in MATLAB 1-based indexing correspond to 4

% Composite Simpson summation
simpson_integral = (dx / 3.0) * sum(weights .* y_vals);

% Reference value resolution
if isempty(exact_reference)
    ref_val = integral(f_handle, a, b, 'RelTol', 1e-12, 'AbsTol', 1e-14);
else
    ref_val = double(exact_reference);
end

truncation_error = abs(simpson_integral - ref_val);

% Assemble solution structure
sol.interval = [a, b];
sol.num_intervals = num_intervals;
sol.step_size_dx = dx;
sol.quadrature_weights = weights;
sol.sampled_nodes_x = x_grid;
sol.sampled_values_y = y_vals;
sol.simpson_value = simpson_integral;
sol.reference_value = ref_val;
sol.absolute_error = truncation_error;
sol.relative_error_pct = (truncation_error / abs(ref_val)) * 100.0;

% Console report
fprintf('=== Simpson''s 1/3 Rule Numerical Quadrature ===\n');
fprintf('Interval Bounds [a, b]      : [%.4f, %.4f]\n', a, b);
fprintf('Subdivision Count (n, even) : %d (\\Delta x = %.6f)\n', num_intervals, dx);
fprintf('Simpson Weight Pattern      : %s\n', mat2str(weights(1:min(end, 9))));
fprintf('Simpson Numerical Integral  : %.10f\n', simpson_integral);
fprintf('Exact / Reference Integral  : %.10f\n', ref_val);
fprintf('Absolute Truncation Error   : %.2e (Rel: %.4f%%)\n', ...
    truncation_error, sol.relative_error_pct);
end