function sol = solve_polar_curve_area(r_sym, theta_bounds)
% SOLVE_POLAR_CURVE_AREA Evaluates the swept enclosed area of a polar function
% r = f(theta) across angular interval [theta1, theta2] using the polar area integral.
% Formula: A = \int_{theta1}^{theta2} \frac{1}{2} * [r(theta)]^2 dtheta
%
% Syntax:
%   sol = solve_polar_curve_area(r_sym, theta_bounds)
%
% Inputs:
%   r_sym        - Symbolic expression r(theta) in variable 'theta' (default: 2*cos(theta))
%   theta_bounds - Two-element angular interval [theta1, theta2] in radians (default: [-pi/2, pi/2])
%
% Outputs:
%   sol          - Struct containing symbolic expression, analytic area, numerical area, and error

arguments
    r_sym        = []
    theta_bounds (1,2) double = [-pi/2, pi/2]
end

syms theta real
if isempty(r_sym), r_sym = 2 * cos(theta); end

th1 = theta_bounds(1);
th2 = theta_bounds(2);

if th1 >= th2
    error('Lower angular limit theta1 must be strictly less than upper limit theta2.');
end

% 1. Formulate Differential Area Integrand: dA = 0.5 * r^2 * dtheta
integrand_sym = simplify(0.5 * (r_sym^2));

% 2. Symbolic Integration
A_analytic_sym = int(integrand_sym, theta, th1, th2);
is_analytic = ~isempty(A_analytic_sym) && ~has(A_analytic_sym, 'int');

if is_analytic
    val_analytic = double(A_analytic_sym);
else
    val_analytic = NaN;
end

% 3. Numerical Validation via Adaptive High-Precision Quadrature
r_handle = matlabFunction(r_sym, 'Vars', theta);
quad_integrand = @(th) 0.5 * (r_handle(th).^2);
val_numerical = integral(quad_integrand, th1, th2, 'RelTol', 1e-12, 'AbsTol', 1e-14);

if is_analytic
    discrepancy = abs(val_analytic - val_numerical);
else
    discrepancy = 0.0;
end

% Assemble solution structure
sol.polar_function = char(r_sym);
sol.area_integrand = char(integrand_sym);
sol.angular_limits = [th1, th2];
sol.angular_span_rad = th2 - th1;
sol.has_analytic_solution = is_analytic;
sol.analytic_area = val_analytic;
sol.numerical_area = val_numerical;
sol.discrepancy = discrepancy;

% Console report
fprintf('=== Polar Coordinate Area Integration Solver ===\n');
fprintf('Radial Function r(\\theta)   : %s\n', sol.polar_function);
fprintf('Area Differential 0.5*r^2  : %s\n', sol.area_integrand);
fprintf('Integration Sector Window  : [%.4f, %.4f] rad (Span: %.2f deg)\n', ...
    th1, th2, rad2deg(sol.angular_span_rad));
if is_analytic
    fprintf('Analytic Enclosed Area (A) : %.10f\n', val_analytic);
else
    fprintf('Analytic Enclosed Area (A) : [Non-Elementary Expression]\n');
end
fprintf('Numerical Quadrature Area  : %.10f (Discrepancy: %.2e)\n', val_numerical, discrepancy);
end