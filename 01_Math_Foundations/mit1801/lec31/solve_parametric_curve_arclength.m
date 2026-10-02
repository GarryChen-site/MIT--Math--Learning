function sol = solve_parametric_curve_arclength(x_sym, y_sym, t_bounds)
% SOLVE_PARAMETRIC_CURVE_ARCLENGTH Computes the exact and numerical arc length
% of a 2D parametric trajectory [x(t), y(t)] over interval [t0, tf] (Lecture 31).
%
% Formulas:
%   Metric Factor : ds/dt = sqrt( (dx/dt)^2 + (dy/dt)^2 )
%   Total Length  : S = \int_t0^tf sqrt( \dot{x}^2 + \dot{y}^2 ) dt
%
% Syntax:
%   sol = solve_parametric_curve_arclength(x_sym, y_sym, t_bounds)
%
% Inputs:
%   x_sym    - Symbolic expression for x(t) in variable 't' (default: t^2)
%   y_sym    - Symbolic expression for y(t) in variable 't' (default: t^3)
%   t_bounds - Two-element parameter vector [t0, tf] (default: [0.0, 1.0])
%
% Outputs:
%   sol      - Struct containing speed metric expression, analytic/numerical length, and cusps

arguments
    x_sym    = []
    y_sym    = []
    t_bounds (1,2) double = [0.0, 1.0]
end

syms t real
if isempty(x_sym), x_sym = t^2; end
if isempty(y_sym), y_sym = t^3; end

t0 = t_bounds(1);
tf = t_bounds(2);

if t0 >= tf
    error('Parameter start bound t0 must be strictly less than end bound tf.');
end

% 1. Derive velocity components
dx_dt = diff(x_sym, t);
dy_dt = diff(y_sym, t);

% Speed metric ds/dt
ds_dt_sym = simplify(sqrt(dx_dt^2 + dy_dt^2));

% Detect potential kinematic cusps/singularities where \dot{x} = 0 and \dot{y} = 0
% Find cusp candidates numerically without invoking solve
cusp_roots = [];
equations = [dx_dt == 0, dy_dt == 0];
initial_guesses = linspace(t0, tf, 25);

for guess = initial_guesses
    try
        root = double(vpasolve(equations, t, guess));
        if isreal(root) && isfinite(root) && root >= t0 && root <= tf
            cusp_roots(end+1) = root; %#ok<AGROW>
        end
    catch
        % Ignore unsuccessful numerical solves
    end
end

valid_cusps = unique(cusp_roots, 'stable');

% 2. Analytic Integration Attempt
S_sym = int(ds_dt_sym, t, t0, tf);
is_analytic = ~isempty(S_sym) && ~has(S_sym, 'int');

if is_analytic
    S_analytic_val = double(S_sym);
else
    S_analytic_val = NaN;
end

% 3. High-Precision Numerical Quadrature
speed_handle = matlabFunction(ds_dt_sym, 'Vars', t);
S_numerical_val = integral(speed_handle, t0, tf, 'RelTol', 1e-12, 'AbsTol', 1e-14);

if is_analytic
    discrepancy = abs(S_analytic_val - S_numerical_val);
else
    discrepancy = 0.0;
end

% Assemble solution structure
sol.x_expression = char(x_sym);
sol.y_expression = char(y_sym);
sol.speed_metric = char(ds_dt_sym);
sol.parameter_bounds = [t0, tf];
sol.cusp_singularities = valid_cusps;
sol.has_analytic_solution = is_analytic;
sol.analytic_arclength = S_analytic_val;
sol.numerical_arclength = S_numerical_val;
sol.discrepancy = discrepancy;

% Console report
fprintf('=== Parametric Curve Arc-Length Solver ===\n');
fprintf('Trajectory Functions       : x(t) = %s | y(t) = %s\n', sol.x_expression, sol.y_expression);
fprintf('Speed Metric ds/dt         : %s\n', sol.speed_metric);
fprintf('Parameter Domain [t0, tf]  : [%.4f, %.4f]\n', t0, tf);
fprintf('Detected Cusp Singularities: %s\n', mat2str(valid_cusps'));
if is_analytic
    fprintf('Analytic Arc Length (S)    : %.10f\n', S_analytic_val);
else
    fprintf('Analytic Arc Length (S)    : [Non-Elementary / Special Function]\n');
end
fprintf('Numerical Quadrature Length: %.10f (Discrepancy: %.2e)\n', S_numerical_val, discrepancy);
end