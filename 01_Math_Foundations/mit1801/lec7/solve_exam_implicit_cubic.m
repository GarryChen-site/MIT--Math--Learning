function sol = solve_exam_implicit_cubic(C_param, x0_val, y0_val)
% SOLVE_EXAM_IMPLICIT_CUBIC Solves implicit derivative and tangent line for
% the cubic algebraic relation: y^3 + 3*x*y^2 = C.
% Formula derived in lecture: dy/dx = -3*y^2 / (3*y^2 + 6*x*y) = -y / (y + 2*x)
%
% Syntax:
%   sol = solve_exam_implicit_cubic(C_param, x0_val, y0_val)
%
% Inputs:
%   C_param - Constant term on RHS (e.g., 8 in Exam 1 review)
%   x0_val  - Point x-coordinate for evaluation (default = 1.0)
%   y0_val  - Point y-coordinate for evaluation (default = 1.0)
%
% Outputs:
%   sol     - Struct containing evaluation status, derivative, and tangent line

arguments
    C_param (1,1) double = 8.0
    x0_val  (1,1) double = 1.0
    y0_val  (1,1) double = 1.0
end

syms x y
F_sym = y^3 + 3*x*y^2 - C_param;

% Partial derivatives
dF_dx = diff(F_sym, x); % 3*y^2
dF_dy = diff(F_sym, y); % 3*y^2 + 6*x*y

% Analytic derivative formula
dy_dx_sym = simplify(-dF_dx / dF_dy);

% Numerical evaluation
dF_dx_val = double(subs(dF_dx, [x, y], [x0_val, y0_val]));
dF_dy_val = double(subs(dF_dy, [x, y], [x0_val, y0_val]));
residual = double(subs(F_sym, [x, y], [x0_val, y0_val]));

if abs(dF_dy_val) < 1e-12
    error('Denominator singularity: dF/dy = 0 leads to infinite tangent slope.');
end

slope_val = -dF_dx_val / dF_dy_val;

% Tangent equation: slope_val * x - y + (y0_val - slope_val * x0_val) = 0
tangent_intercept = y0_val - slope_val * x0_val;

sol.C = C_param;
sol.x0 = x0_val;
sol.y0 = y0_val;
sol.equation_residual = residual;
sol.symbolic_derivative = char(dy_dx_sym);
sol.tangent_slope = slope_val;
sol.tangent_y_intercept = tangent_intercept;

% Console output
fprintf('=== Exam 1 Implicit Cubic Solver: y^3 + 3xy^2 = %.2f ===\n', C_param);
fprintf('Operating Point (x0, y0)   : (%.4f, %.4f) [Residual: %.2e]\n', x0_val, y0_val, residual);
fprintf('Symbolic Formula (dy/dx)   : %s\n', sol.symbolic_derivative);
fprintf('Evaluated Slope (dy/dx)    : %.6f\n', slope_val);
fprintf('Tangent Line Formulation   : y = %.4f * x + (%.4f)\n', slope_val, tangent_intercept);
end