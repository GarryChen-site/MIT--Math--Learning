function sol = solve_golf_ball_surface_area(R_ball, r_dimple, N_dimples, pi_approx)
% SOLVE_GOLF_BALL_SURFACE_AREA Solves the surface area of a golf ball with hemispherical
% dimples using three levels of spherical cap approximations (PS2B Part II Problem 1).
% Models:
%   i)   Flat disk removed   : A = 4*pi*R^2 + N * pi * r^2
%   ii)  2nd-order expansion : A = 4*pi*R^2 + N * pi * r^2 - N * (pi * r^4) / (4 * R^2)
%   iii) Exact spherical cap : A = 4*pi*R^2 - N * 2*pi*R^2*(1 - sqrt(1 - (r/R)^2)) + N * 2*pi*r^2
%
% Syntax:
%   sol = solve_golf_ball_surface_area(R_ball, r_dimple, N_dimples, pi_approx)
%
% Inputs:
%   R_ball    - Sphere radius in cm (default = 1.5 cm, diameter 3 cm)
%   r_dimple  - Dimple radius in cm (default = 0.15 cm, diameter 3 mm)
%   N_dimples - Number of dimples (default = 100)
%   pi_approx - Value of pi to use (default = 3.14 as specified in contest)
%
% Outputs:
%   sol       - Struct containing computed areas, error discrepancies, and packing limit

arguments
    R_ball    (1,1) double {mustBePositive} = 1.5
    r_dimple  (1,1) double {mustBePositive} = 0.15
    N_dimples (1,1) double {mustBeInteger, mustBePositive} = 100
    pi_approx (1,1) double {mustBePositive} = 3.14
end

if r_dimple >= R_ball
    error('Dimple radius r must be strictly less than golf ball radius R.');
end

ratio_u = r_dimple / R_ball;

% 1. Model (i): Flat Disk Removed Approximation
% Area removed per dimple approx pi*r^2. Added hemispherical area = 2*pi*r^2.
% Net change per dimple = + pi*r^2
area_base_sphere = 4 * pi_approx * (R_ball^2);
area_flat_approx = area_base_sphere + N_dimples * (pi_approx * (r_dimple^2));

% 2. Model (ii): 2nd-Order Curvature Taylor Approximation
% Cap area approx pi*r^2 + (pi*r^4)/(4*R^2)
% Net = 4*pi*R^2 + N*pi*r^2 - N*(pi*r^4)/(4*R^2)
curvature_correction = N_dimples * (pi_approx * (r_dimple^4)) / (4 * (R_ball^2));
area_quad_approx = area_flat_approx - curvature_correction;

% 3. Model (iii): Exact Analytic Closed-Form
cap_area_exact = 2 * pi_approx * (R_ball^2) * (1 - sqrt(1 - ratio_u^2));
hemisphere_added = 2 * pi_approx * (r_dimple^2);
area_exact = area_base_sphere - N_dimples * cap_area_exact + N_dimples * hemisphere_added;

% 4. Extra credit: Maximum non-overlapping dimples packing estimate (Part d)
% Theoretical max hexagonal circle packing density on a sphere is ~ 0.9069
% Cap surface area fraction per dimple:
packing_fraction = 0.9069;
max_dimples_theoretical = floor((area_base_sphere * packing_fraction) / cap_area_exact);

% Assemble solution structure
sol.R_ball = R_ball;
sol.r_dimple = r_dimple;
sol.N_dimples = N_dimples;
sol.ratio_r_over_R = ratio_u;
sol.area_flat_approx = area_flat_approx;
sol.area_flat_rounded = round(area_flat_approx, 2);
sol.area_quad_approx = area_quad_approx;
sol.area_quad_rounded = round(area_quad_approx, 2);
sol.area_exact = area_exact;
sol.area_exact_rounded = round(area_exact, 2);
sol.flat_discrepancy = area_flat_approx - area_exact;
sol.quad_discrepancy = area_quad_approx - area_exact;
sol.max_dimples_estimate = max_dimples_theoretical;

% Console report
fprintf('=== PS2B Problem 1: Golf Ball Surface Area Analysis ===\n');
fprintf('Geometry: R = %.2f cm, r = %.3f cm (Ratio r/R = %.3f), N = %d\n', ...
    R_ball, r_dimple, ratio_u, N_dimples);
fprintf('Model (i)   [Flat Disk Approximation] : %.5f cm^2 -> Rounded: %.2f cm^2\n', ...
    area_flat_approx, sol.area_flat_rounded);
fprintf('Model (ii)  [2nd-Order Approximated]  : %.5f cm^2 -> Rounded: %.2f cm^2\n', ...
    area_quad_approx, sol.area_quad_rounded);
fprintf('Model (iii) [Exact Closed-Form]       : %.5f cm^2 -> Rounded: %.2f cm^2\n', ...
    area_exact, sol.area_exact_rounded);
fprintf('Contest Result Verification: Flat Disk gives %.2f vs Exact %.2f (FLAT DISK FAILS!)\n', ...
    sol.area_flat_rounded, sol.area_exact_rounded);
fprintf('Extra Credit: Theoretical Max Dimple Capacity : ~%d dimples\n', max_dimples_theoretical);
end