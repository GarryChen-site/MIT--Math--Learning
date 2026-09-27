function sol = solve_enclosed_area_horizontal(c_parabola, m_line, b_line)
% SOLVE_ENCLOSED_AREA_HORIZONTAL Computes the exact and numerical enclosed area
% between a sideways parabola x = c*y^2 and straight line x = m*y + b.
% Integration: A = \int_{y1}^{y2} [ (m*y + b) - (c*y^2) ] dy
%
% Syntax:
%   sol = solve_enclosed_area_horizontal(c_parabola, m_line, b_line)
%
% Inputs:
%   c_parabola - Quadratic coefficient for x = c*y^2 (default = 1.0)
%   m_line     - Slope parameter for x = m*y + b (default = 1.0)
%   b_line     - Intercept parameter for x = m*y + b (default = 2.0)
%
% Outputs:
%   sol        - Struct containing intersection points, integration limits, and areas

arguments
    c_parabola (1,1) double {mustBePositive} = 1.0
    m_line     (1,1) double = 1.0
    b_line     (1,1) double = 2.0
end

% Intersection equation: c*y^2 - m*y - b = 0
discriminant = m_line^2 + 4 * c_parabola * b_line;
if discriminant <= 0
    error('Curves do not intersect or intersect at a single point; no enclosed area exists.');
end

% 1. Intersections in y-domain (lower limit y1, upper limit y2)
y1 = (m_line - sqrt(discriminant)) / (2 * c_parabola);
y2 = (m_line + sqrt(discriminant)) / (2 * c_parabola);

% Associated x-coordinates
x1 = c_parabola * (y1^2);
x2 = c_parabola * (y2^2);

% 2. Exact Integration via horizontal slicing (right curve minus left curve):
% Integrand: g(y) = (m*y + b) - c*y^2
% Antiderivative: G(y) = m*y^2 / 2 + b*y - c*y^3 / 3
G_anti = @(y) (m_line * y.^2) / 2.0 + b_line * y - (c_parabola * y.^3) / 3.0;
area_analytic = G_anti(y2) - G_anti(y1);

% 3. Numerical Verification via adaptive quadrature
integrand_h = @(y) (m_line * y + b_line) - c_parabola * (y.^2);
area_numerical = integral(integrand_h, y1, y2, 'RelTol', 1e-12, 'AbsTol', 1e-14);
discrepancy = abs(area_analytic - area_numerical);

% Assemble solution structure
sol.c_parabola = c_parabola;
sol.m_line = m_line;
sol.b_line = b_line;
sol.intersections_y = [y1, y2];
sol.intersections_cartesian = [x1, y1; x2, y2];
sol.analytic_area = area_analytic;
sol.numerical_area = area_numerical;
sol.discrepancy = discrepancy;

% Console report
fprintf('=== Enclosed Area Solver via Horizontal Slices ===\n');
fprintf('Curves: Parabola x = %.2f*y^2 | Line x = %.2f*y + %.2f\n', c_parabola, m_line, b_line);
fprintf('Intersection 1 (Lower)    : (x = %.4f, y = %.4f)\n', x1, y1);
fprintf('Intersection 2 (Upper)    : (x = %.4f, y = %.4f)\n', x2, y2);
fprintf('Analytic Enclosed Area    : %.10f\n', area_analytic);
fprintf('Numerical Quadrature Area : %.10f\n', area_numerical);
fprintf('Integration Discrepancy   : %.2e\n', discrepancy);
end