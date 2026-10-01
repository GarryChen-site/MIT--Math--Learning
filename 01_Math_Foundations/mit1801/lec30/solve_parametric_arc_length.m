function sol = solve_parametric_arc_length(x_sym, y_sym, t_bounds, num_resample_pts)
% SOLVE_PARAMETRIC_ARC_LENGTH Computes the exact and numerical arc length
% for a 2D parametric trajectory [x(t), y(t)] and generates uniformly spaced
% arc-length waypoints for Frenet-frame path tracking (Lecture 30 Example 5-7).
%
% Formulas:
%   Speed / Metric Factor : ds/dt = sqrt( (dx/dt)^2 + (dy/dt)^2 )
%   Total Arc Length      : S = \int_{t0}^{tf} sqrt( \dot{x}^2 + \dot{y}^2 ) dt
%
% Syntax:
%   sol = solve_parametric_arc_length(x_sym, y_sym, t_bounds, num_resample_pts)
%
% Inputs:
%   x_sym            - Symbolic expression for x(t) in variable 't' (default: 2*sin(t))
%   y_sym            - Symbolic expression for y(t) in variable 't' (default: cos(t))
%   t_bounds         - Two-element parameter vector [t0, tf] (default: [0.0, 2*pi])
%   num_resample_pts - Number of equidistant arc-length waypoints (default = 50)
%
% Outputs:
%   sol              - Struct containing total length, speed function, and resampled nodes

arguments
    x_sym            = []
    y_sym            = []
    t_bounds         (1,2) double = [0.0, 2*pi]
    num_resample_pts (1,1) double {mustBeInteger, mustBePositive} = 50
end

syms t real
if isempty(x_sym), x_sym = 2 * sin(t); end
if isempty(y_sym), y_sym = cos(t); end

t0 = t_bounds(1);
tf = t_bounds(2);

if t0 >= tf
    error('Lower bound t0 must be strictly less than upper bound tf.');
end

% 1. Symbolic derivatives for tangent velocities
dx_dt_sym = diff(x_sym, t);
dy_dt_sym = diff(y_sym, t);

% Speed metric ds/dt
speed_sym = simplify(sqrt(dx_dt_sym^2 + dy_dt_sym^2));

% 2. Numerical Integration via Adaptive High-Precision Quadrature
speed_handle = matlabFunction(speed_sym, 'Vars', t);
total_length = integral(speed_handle, t0, tf, 'RelTol', 1e-12, 'AbsTol', 1e-14);

% 3. Uniform Arc-Length Resampling for Waypoint Tracking
t_dense = linspace(t0, tf, 2000);
speed_dense = speed_handle(t_dense);
s_dense = cumtrapz(t_dense, speed_dense);

% Target equidistant arc lengths
s_uniform = linspace(0, s_dense(end), num_resample_pts);
% Invert s(t) => t(s) via shape-preserving monotonic interpolation
t_uniform = interp1(s_dense, t_dense, s_uniform, 'pchip');

x_func = matlabFunction(x_sym, 'Vars', t);
y_func = matlabFunction(y_sym, 'Vars', t);

x_waypoints = x_func(t_uniform);
y_waypoints = y_func(t_uniform);

% Assemble solution structure
sol.x_expression = char(x_sym);
sol.y_expression = char(y_sym);
sol.parameter_bounds = [t0, tf];
sol.speed_profile_expr = char(speed_sym);
sol.total_arc_length = total_length;
sol.resampled_arc_lengths = s_uniform';
sol.resampled_parameter_t = t_uniform';
sol.resampled_waypoints = [x_waypoints', y_waypoints'];

% Console report
fprintf('=== Parametric Curve Arc-Length & Path Resampling Solver ===\n');
fprintf('Parametric Function x(t)    : %s\n', sol.x_expression);
fprintf('Parametric Function y(t)    : %s\n', sol.y_expression);
fprintf('Parameter Horizon [t0, tf]  : [%.4f, %.4f]\n', t0, tf);
fprintf('Speed Metric Factor ds/dt   : %s\n', sol.speed_profile_expr);
fprintf('Total Accumulated Arc Length: %.10f units\n', total_length);
fprintf('Generated Equidistant Nodes : %d Waypoints (Constant \\Delta s = %.6f)\n', ...
    num_resample_pts, total_length / (num_resample_pts - 1));
end