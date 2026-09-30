function report = solve_circle_segment_trig_sub(radius_a, upper_limit_x)
% SOLVE_CIRCLE_SEGMENT_TRIG_SUB Evaluates the circular segment integral:
% A(x) = \int_0^x sqrt(a^2 - t^2) dt via trigonometric substitution t = a*sin(u).
% Formula: A(x) = (a^2 / 2) * asin(x / a) + (x / 2) * sqrt(a^2 - x^2)
% Geometric breakdown: Circular Sector Area + Right Triangle Area.
%
% Syntax:
%   report = solve_circle_segment_trig_sub(radius_a, upper_limit_x)
%
% Inputs:
%   radius_a      - Circle radius a > 0 (default = 2.0)
%   upper_limit_x - Integration bound x with 0 <= x <= a (default = 1.0)
%
% Outputs:
%   report        - Struct containing geometric parts, total area, and numerical validation

arguments
    radius_a      (1,1) double {mustBePositive} = 2.0
    upper_limit_x (1,1) double {mustBeNonnegative} = 1.0
end

if upper_limit_x > radius_a
    error('Upper limit x (%.4f) cannot exceed circle radius a (%.4f).', upper_limit_x, radius_a);
end

% 1. Substitution Parameters: t = a*sin(u) => u = asin(x / a)
u_limit = asin(upper_limit_x / radius_a);
u_limit_deg = rad2deg(u_limit);

% 2. Geometric Component Breakdown:
% Part A: Circular sector area corresponding to angle u: 0.5 * a^2 * u
sector_area = 0.5 * (radius_a^2) * u_limit;

% Part B: Right triangle area under chord: 0.5 * base * height = 0.5 * x * sqrt(a^2 - x^2)
triangle_height = sqrt(radius_a^2 - upper_limit_x^2);
triangle_area = 0.5 * upper_limit_x * triangle_height;

% Total Analytic Area via Trig Substitution:
area_analytic = sector_area + triangle_area;

% 3. Numerical Integration Verification
integrand_circle = @(t) sqrt(radius_a^2 - t.^2);
area_numerical = integral(integrand_circle, 0, upper_limit_x, 'RelTol', 1e-12, 'AbsTol', 1e-14);

discrepancy = abs(area_analytic - area_numerical);

% Assemble structured report
report.radius = radius_a;
report.upper_limit_x = upper_limit_x;
report.substitution_angle_rad = u_limit;
report.substitution_angle_deg = u_limit_deg;
report.circular_sector_area = sector_area;
report.right_triangle_area = triangle_area;
report.total_analytic_area = area_analytic;
report.numerical_area_quadrature = area_numerical;
report.discrepancy = discrepancy;

% Console report
fprintf('=== Circular Segment Area via Trigonometric Substitution ===\n');
fprintf('Circle Radius (a)           : %.4f | Integration Limit (x): %.4f\n', radius_a, upper_limit_x);
fprintf('Substituted Angle u=asin(x/a): %.4f rad (%.2f deg)\n', u_limit, u_limit_deg);
fprintf('Circular Sector Area (A1)   : %.8f\n', sector_area);
fprintf('Right Triangle Area (A2)    : %.8f\n', triangle_area);
fprintf('Total Analytic Area (A1+A2) : %.10f\n', area_analytic);
fprintf('Numerical Integral Check    : %.10f (Discrepancy: %.2e)\n', area_numerical, discrepancy);
end