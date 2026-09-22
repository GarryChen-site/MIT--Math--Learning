function sol = solve_implicit_derivative(F_sym, x0_val, y0_val)
% SOLVE_IMPLICIT_DERIVATIVE Computes the symbolic and numerical derivative
% dy/dx for an arbitrary implicit equation F(x, y) = 0 using partial differentials.
% Formula: dy/dx = - (dF/dx) / (dF/dy)
%
% Syntax:
%   sol = solve_implicit_derivative(F_sym, x0_val, y0_val)
%
% Inputs:
%   F_sym  - Symbolic expression representing F(x, y) = 0 (using variables x, y)
%   x0_val - Operating evaluation coordinate x0
%   y0_val - Operating evaluation coordinate y0
%
% Outputs:
%   sol    - Struct containing symbolic partials, dy/dx formula, evaluated slope,
%            and tangent line coefficients [A, B, C] for A*x + B*y + C = 0

arguments
    F_sym  = []
    x0_val (1,1) double = 0.6
    y0_val (1,1) double = 0.8
end

syms x y
if isempty(F_sym), F_sym = x^2 + y^2 - 1; end

% 1. Partial differentiation
dF_dx = diff(F_sym, x);
dF_dy = diff(F_sym, y);

% 2. Implicit derivative formula
dy_dx_sym = -dF_dx / dF_dy;

% 3. Check equation satisfaction at given point
residual = double(subs(F_sym, [x, y], [x0_val, y0_val]));
if abs(residual) > 1e-4
    warning('The point (%.4f, %.4f) does not strictly lie on curve F(x,y)=0 (Residual: %.2e)', ...
        x0_val, y0_val, residual);
end

% 4. Numerical evaluation
dF_dx_num = double(subs(dF_dx, [x, y], [x0_val, y0_val]));
dF_dy_num = double(subs(dF_dy, [x, y], [x0_val, y0_val]));

if abs(dF_dy_num) < 1e-12
    error('Singularity detected: dF/dy is zero (vertical tangent / kinematic singularity).');
end

slope_val = -dF_dx_num / dF_dy_num;

% Tangent line: y - y0 = m * (x - x0) => m*x - y + (y0 - m*x0) = 0
tangent_line_coeffs = [slope_val, -1, (y0_val - slope_val * x0_val)];

% Assemble solution structure
sol.x0 = x0_val;
sol.y0 = y0_val;
sol.residual = residual;
sol.symbolic_F = char(F_sym);
sol.symbolic_dy_dx = char(simplify(dy_dx_sym));
sol.partial_dF_dx = dF_dx_num;
sol.partial_dF_dy = dF_dy_num;
sol.tangent_slope = slope_val;
sol.tangent_line_coeffs = tangent_line_coeffs;

% Console output
fprintf('=== Implicit Differentiation Analysis ===\n');
fprintf('Implicit Equation F(x, y) = 0 : %s = 0\n', sol.symbolic_F);
fprintf('Evaluated Operating Point     : (x0 = %.4f, y0 = %.4f)\n', x0_val, y0_val);
fprintf('Symbolic Derivative dy/dx     : %s\n', sol.symbolic_dy_dx);
fprintf('Partial Derivatives at Point  : dF/dx = %.4f | dF/dy = %.4f\n', dF_dx_num, dF_dy_num);
fprintf('Tangent Line Slope (dy/dx)    : %.6f\n', slope_val);
fprintf('Tangent Line Equation         : %.4f*x + %.4f*y + %.4f = 0\n', ...
    tangent_line_coeffs(1), tangent_line_coeffs(2), tangent_line_coeffs(3));
end