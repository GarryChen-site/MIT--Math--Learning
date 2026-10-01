function sol = solve_circle_segment_integral(radius_a, upper_bound_x)
% SOLVE_CIRCLE_SEGMENT_INTEGRAL Computes the exact area under a circular arc
% via inverse trigonometric substitution t = a*sin(u) (Lecture 28).
%
% Formulas:
%   Integral : \int_0^x sqrt(a^2 - t^2) dt
%   Analytic : (a^2 / 2) * asin(x / a) + (x / 2) * sqrt(a^2 - x^2)
%   Breakdown: Sector Area [ (a^2 / 2) * u ] + Triangle Area [ (x / 2) * y ]
%
% Syntax:
%   sol = solve_circle_segment_integral(radius_a, upper_bound_x)
%
% Inputs:
%   radius_a      - Circle radius a > 0 (default = 2.0)
%   upper_bound_x - Upper boundary x satisfying 0 <= x <= a (default = 1.0)
%
% Outputs:
%   sol           - Struct containing sector area, triangle area, total area, and error

arguments
    radius_a      (1,1) double {mustBePositive} = 2.0
    upper_bound_x (1,1) double {mustBeNonnegative} = 1.0
end

if upper_bound_x > radius_a
    error('Upper limit x (%.4f) must not exceed circle radius a (%.4f).', upper_bound_x, radius_a);
end

% 1. Inverse Substitution Angle: u = asin(x / a)
u_angle = asin(upper_bound_x / radius_a);

% 2. Geometric Component Decomposition:
% Part 1: Circular Sector Area
sector_val = 0.5 * (radius_a^2) * u_angle;

% Part 2: Right Triangle Area (base = x, height = sqrt(a^2 - x^2))
tri_height = sqrt(radius_a^2 - upper_bound_x^2);
triangle_val = 0.5 * upper_bound_x * tri_height;

% Total Analytic Area
total_analytic = sector_val + triangle_val;

% 3. Numerical Validation via Adaptive Quadrature
f_circ = @(t) sqrt(radius_a^2 - t.^2);
total_numerical = integral(f_circ, 0, upper_bound_x, 'RelTol', 1e-12, 'AbsTol', 1e-14);

discrepancy = abs(total_analytic - total_numerical);

% Assemble solution structure
sol.circle_radius = radius_a;
sol.upper_limit = upper_bound_x;
sol.substitution_angle_rad = u_angle;
sol.substitution_angle_deg = rad2deg(u_angle);
sol.sector_area_component = sector_val;
sol.triangle_area_component = triangle_val;
sol.total_analytic_area = total_analytic;
sol.total_numerical_area = total_numerical;
sol.discrepancy = discrepancy;

% Console report
fprintf('=== Circular Arc Area via Inverse Substitution ===\n');
fprintf('Circle Radius (a)         : %.4f | Integration Limit (x): %.4f\n', radius_a, upper_bound_x);
fprintf('Inverse Angle u = asin(x/a): %.6f rad (%.2f deg)\n', u_angle, sol.substitution_angle_deg);
fprintf('Sector Component (0.5 a^2 u): %.8f\n', sector_val);
fprintf('Triangle Component (0.5 x y): %.8f (Height y = %.4f)\n', triangle_val, tri_height);
fprintf('Total Analytic Area         : %.10f\n', total_analytic);
fprintf('Numerical Quadrature Check  : %.10f (Error: %.2e)\n', total_numerical, discrepancy);
end