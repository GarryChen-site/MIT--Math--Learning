function report = solve_wire_cut_extreme_area(total_wire_length, num_sides)
% SOLVE_WIRE_CUT_EXTREME_AREA Solves optimal cut allocation of a wire of length L
% into two regular polygons with N sides.
% Formulates total area: A(x) = C_shape * [ (x/N)^2 + ((L - x)/N)^2 ]
% Evaluates critical points, second-derivative curvature, and boundary limits x in [0, L].
%
% Syntax:
%   report = solve_wire_cut_extreme_area(total_wire_length, num_sides)
%
% Inputs:
%   total_wire_length - Total length of wire L (default = 1.0)
%   num_sides         - Number of polygon sides N (default = 4 for squares)
%
% Outputs:
%   report            - Struct containing critical and boundary coordinates and global classifications

arguments
    total_wire_length (1,1) double {mustBePositive} = 1.0
    num_sides         (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(num_sides, 3)} = 4
end

% Area constant for regular N-gon with side length s:
% A_polygon = (N / 4) * cot(pi / N) * s^2
k_shape = (num_sides / 4.0) * cot(pi / num_sides);

% Segment perimeter lengths are x and (L - x)
% Side lengths are s1 = x/N, s2 = (L - x)/N
% Total Area A(x) = (k_shape / N^2) * [ x^2 + (L - x)^2 ]
coeff_total = k_shape / (num_sides^2);

% 1. Critical Point Analysis:
% dA/dx = coeff_total * [ 2*x - 2*(L - x) ] = 0 => 4*x - 2*L = 0 => x* = L / 2
x_crit = total_wire_length / 2.0;
A_crit = coeff_total * (x_crit^2 + (total_wire_length - x_crit)^2);


% 2. Boundary Endpoints Analysis (x = 0 and x = L):
% All wire dedicated to a single regular polygon
A_boundary = coeff_total * (0^2 + total_wire_length^2);

% Classify global extrema
global_min_area = A_crit;
global_min_location = x_crit;

global_max_area = A_boundary;
global_max_locations = [0.0, total_wire_length];

% Assemble structured report
report.total_wire_length = total_wire_length;
report.polygon_sides = num_sides;
report.critical_cut_x = x_crit;
report.critical_area = A_crit;
report.critical_point_type = 'Strict Local Minimum (d2A/dx2 > 0)';
report.boundary_area = A_boundary;
report.global_minimum = global_min_area;
report.global_min_cut_point = global_min_location;
report.global_maximum = global_max_area;
report.global_max_cut_points = global_max_locations;

% Console report
fprintf('=== Wire Partitioning & Enclosed Area Optimization ===\n');
fprintf('Total Wire Length (L)      : %.4f\n', total_wire_length);
fprintf('Polygon Configuration      : Two Regular %d-gons\n', num_sides);
fprintf('Critical Cut Position (x*) : %.4f (Equal split: L/2)\n', x_crit);
fprintf('Critical Point Area        : %.6f [%s]\n', A_crit, report.critical_point_type);
fprintf('Boundary Enclosed Area     : %.6f (Single undivided polygon)\n', A_boundary);
fprintf('GLOBAL MINIMUM             : Area = %.6f (Occurs at symmetric center x = %.4f)\n', ...
    global_min_area, global_min_location);
fprintf('GLOBAL MAXIMUM             : Area = %.6f (Occurs at saturation bounds x = 0 or x = %.4f)\n', ...
    global_max_area, total_wire_length);
end