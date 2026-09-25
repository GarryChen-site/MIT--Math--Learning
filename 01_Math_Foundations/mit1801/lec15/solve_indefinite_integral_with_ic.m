function result = solve_indefinite_integral_with_ic(f_sym, x0_val, y0_val)
% SOLVE_INDEFINITE_INTEGRAL_WITH_IC Solves indefinite integral F(x) = int f(x)dx,
% determines unique integration constant c matching initial condition F(x0) = y0,
% and analytically validates derivative reconstruction F'(x) == f(x).
%
% Syntax:
%   result = solve_indefinite_integral_with_ic(f_sym, x0_val, y0_val)
%
% Inputs:
%   f_sym  - Symbolic integrand in variable 'x' (default: x^3 * (x^4 + 2)^5)
%   x0_val - Initial state independent variable coordinate (default: 0.0)
%   y0_val - Initial state accumulated condition (default: 1.0)
%
% Outputs:
%   result - Struct containing general antiderivative, calibrated constant c, and exact solution

arguments
    f_sym  = []
    x0_val (1,1) double = 0.0
    y0_val (1,1) double = 1.0
end

syms x c_const real
if isempty(f_sym), f_sym = x^3 * (x^4 + 2)^5; end

% 1. Indefinite integration
F_general_sym = int(f_sym, x);

% 2. Solve for integration constant matching initial condition F(x0) + c = y0
F_at_x0 = subs(F_general_sym, x, x0_val);
c_solved = solve(F_at_x0 + c_const == y0_val, c_const);
c_val = double(c_solved);

% Particular solution satisfying initial state
F_particular_sym = F_general_sym + c_val;

% 3. Analytic verification: d/dx [F_particular] == f(x)
dF_dx_sym = diff(F_particular_sym, x);
verification_residual = simplify(dF_dx_sym - f_sym);
is_derivative_verified = isequal(verification_residual, sym(0));

% Assemble solution structure
result.integrand = char(f_sym);
result.general_antiderivative = char(F_general_sym);
result.initial_condition_point = [x0_val, y0_val];
result.resolved_constant_c = c_val;
result.particular_solution = char(F_particular_sym);
result.is_derivative_verified = is_derivative_verified;

% Console report
fprintf('=== Antiderivative Solver with Initial Condition ===\n');
fprintf('Integrand f(x)              : %s\n', result.integrand);
fprintf('General Antiderivative F(x) : %s + c\n', result.general_antiderivative);
fprintf('Initial Condition [x0, y0]  : [%.4f, %.4f]\n', x0_val, y0_val);
fprintf('Calibrated Constant c       : %.6f\n', c_val);
fprintf('Particular State Solution   : %s\n', result.particular_solution);
fprintf('Derivative Check (F'' == f)  : %s\n', mat2str(is_derivative_verified));
end