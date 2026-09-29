function report = solve_revolution_volume_decision(f_sym, x_bounds, axis_of_revolution)
% SOLVE_REVOLUTION_VOLUME_DECISION Implements Exam 3 Review decision rubric for
% solids of revolution: evaluates volumes via both Disks/Washers and Cylindrical Shells.
%
% Options for axis_of_revolution: 'x' or 'y'
% If rotated around y-axis:
%   Shells : V = 2*pi * \int x * f(x) dx
%   Disks  : V = pi * \int (x_outer^2 - x_inner^2) dy
%
% Syntax:
%   report = solve_revolution_volume_decision(f_sym, x_bounds, axis_of_revolution)
%
% Inputs:
%   f_sym              - Upper bound function y = f(x) in variable 'x' (default: x^2)
%   x_bounds           - Radial limits [x1, x2] with 0 <= x1 < x2 (default: [0.0, 1.0])
%   axis_of_revolution - Axis char: 'x' or 'y' (default: 'y' for witch's cauldron)
%
% Outputs:
%   report             - Struct containing chosen primary method, dual volume checks, and error

arguments
    f_sym              = []
    x_bounds           (1,2) double = [0.0, 1.0]
    axis_of_revolution (1,1) char {mustBeMember(axis_of_revolution, ['x', 'y'])} = 'y'
end

syms x y real
if isempty(f_sym), f_sym = x^2; end

x1 = x_bounds(1);
x2 = x_bounds(2);

if x1 < 0 || x1 >= x2
    error('Invalid interval: Require 0 <= x1 < x2.');
end

y1 = double(subs(f_sym, x, x1));
y2 = double(subs(f_sym, x, x2));

if axis_of_revolution == 'y'
    % Rotation around Y-Axis
    primary_method = 'Cylindrical Shell Method (Integrate dx along radius x)';

    % Method A: Shells -> V = 2*pi * \int_{x1}^{x2} x * f(x) dx
    integrand_shell = simplify(2 * pi * x * f_sym);
    V_shell_sym = int(integrand_shell, x, x1, x2);
    V_shell = double(V_shell_sym);

    % Method B: Disks/Washers -> V = pi * \int_{y1}^{y2} (x_outer^2 - x_inner^2) dy
    % Invert y = f(x) => x = g(y)
    x_inv = solve(f_sym == y, x);
    x_inv_positive = x_inv(x_inv >= 0); % Select positive branch

    % Disk volume of outer cylinder minus inner paraboloid
    % V_disk = pi * x2^2 * (y2 - y1) - pi * \int_{y1}^{y2} [g(y)]^2 dy
    V_outer_cylinder = pi * (x2^2) * (y2 - y1);
    integrand_disk = simplify(pi * (x_inv_positive(1)^2));
    V_inner_bowl = double(int(integrand_disk, y, y1, y2));
    V_disk = V_outer_cylinder - V_inner_bowl;

    dual_discrepancy = abs(V_shell - V_disk);

else
    % Rotation around X-Axis
    primary_method = 'Disk Method (Integrate dx along cross-section \pi*y^2)';

    % Method A: Disks -> V = pi * \int_{x1}^{x2} [f(x)]^2 dx
    integrand_disk = simplify(pi * (f_sym^2));
    V_disk_sym = int(integrand_disk, x, x1, x2);
    V_disk = double(V_disk_sym);

    % Method B: Shells -> V = 2*pi * \int_{y1}^{y2} y * (x2 - g(y)) dy
    x_inv = solve(f_sym == y, x);
    x_inv_positive = x_inv(x_inv >= 0);
    integrand_shell = simplify(2 * pi * y * (x2 - x_inv_positive(1)));
    V_shell = double(int(integrand_shell, y, y1, y2));

    dual_discrepancy = abs(V_disk - V_shell);
end

% Assemble structured report
report.profile_function = char(f_sym);
report.rotation_axis = axis_of_revolution;
report.recommended_method = primary_method;
report.x_bounds = [x1, x2];
report.y_bounds = [y1, y2];
report.volume_by_shells = V_shell;
report.volume_by_disks = V_disk;
report.dual_path_discrepancy = dual_discrepancy;
report.is_dual_verified = (dual_discrepancy < 1e-8);

% Console report
fprintf('=== Solids of Revolution: Algorithmic Decision Solver ===\n');
fprintf('Profile Function y = f(x)  : %s\n', report.profile_function);
fprintf('Rotation Axis              : %c-Axis\n', axis_of_revolution);
fprintf('Recommended Slicing Tool   : %s\n', primary_method);
fprintf('Domain Limits [x1, x2]     : [%.4f, %.4f] => y in [%.4f, %.4f]\n', x1, x2, y1, y2);
fprintf('Volume via Shells Method   : %.10f\n', V_shell);
fprintf('Volume via Disks Method    : %.10f\n', V_disk);
fprintf('Dual Methods Consistency   : %s (Discrepancy: %.2e)\n', ...
    mat2str(report.is_dual_verified), dual_discrepancy);
end