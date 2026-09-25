function sol = solve_separable_first_order_ode(f_sym, g_sym, x0_val, y0_val, x_target)
% SOLVE_SEPARABLE_FIRST_ORDER_ODE Solves initial value problem (IVP) for separable ODE:
% dy/dx = f(x) * g(y) with condition y(x0) = y0.
% Analytical formulation: \int (1/g(y)) dy = \int f(x) dx + c
% Compares analytic closed-form solution with numerical Runge-Kutta integration (ode45).
%
% Syntax:
%   sol = solve_separable_first_order_ode(f_sym, g_sym, x0_val, y0_val, x_target)
%
% Inputs:
%   f_sym    - Symbolic term f(x) depending only on 'x' (default: -x)
%   g_sym    - Symbolic term g(y) depending only on 'y' (default: y)
%   x0_val   - Initial state independent coordinate (default = 0.0)
%   y0_val   - Initial state dependent value (default = 1.0)
%   x_target - Evaluation destination coordinate for accuracy check (default = 1.5)
%
% Outputs:
%   sol      - Struct containing symbolic solutions, integration constant c, and error

arguments
    f_sym    = []
    g_sym    = []
    x0_val   (1,1) double = 0.0
    y0_val   (1,1) double = 1.0
    x_target (1,1) double = 1.5
end

syms x y c_const real
if isempty(f_sym), f_sym = -x; end
if isempty(g_sym), g_sym = y; end

% 1. Indefinite integration on separated sides
% Left side: H(y) = \int (1 / g(y)) dy
% Right side: F(x) = \int f(x) dx
H_sym = int(1 / g_sym, y);
F_sym = int(f_sym, x);

% 2. Solve for unique constant c matching initial state [x0, y0]
H_at_y0 = subs(H_sym, y, y0_val);
F_at_x0 = subs(F_sym, x, x0_val);
c_solved = double(solve(H_at_y0 == F_at_x0 + c_const, c_const));

% Implicit equation: H(y) = F(x) + c
implicit_eqn = (H_sym == F_sym + c_solved);

% Attempt explicit symbolic inversion y = H^(-1)(F(x) + c)
y_explicit_sym = solve(implicit_eqn, y);
if isempty(y_explicit_sym)
    y_analytic_func = @(xv) NaN;
    y_analytic_target = NaN;
else
    % Select branch matching initial condition
    y_branch = y_explicit_sym(1);
    y_analytic_target = double(subs(y_branch, x, x_target));
    y_analytic_func = matlabFunction(y_branch);
end

% 3. Numerical ODE integration via ode45
ode_rhs = @(x_val, y_val) double(subs(f_sym * g_sym, [x, y], [x_val, y_val]));
t_span = [x0_val, x_target];
[~, Y_num] = ode45(ode_rhs, t_span, y0_val);
y_num_target = Y_num(end);

discrepancy = abs(y_analytic_target - y_num_target);

% Assemble solution structure
sol.x0 = x0_val;
sol.y0 = y0_val;
sol.target_x = x_target;
sol.separated_H_y = char(H_sym);
sol.separated_F_x = char(F_sym);
sol.integration_constant_c = c_solved;
sol.analytic_solution = char(y_explicit_sym(1));
sol.analytic_target_value = y_analytic_target;
sol.numerical_target_value = y_num_target;
sol.absolute_discrepancy = discrepancy;

% Console report
fprintf('=== Separable ODE Analytic & Numerical Solver ===\n');
fprintf('ODE Formulation            : dy/dx = (%s) * (%s)\n', char(f_sym), char(g_sym));
fprintf('Initial Condition [x0, y0] : [%.4f, %.4f]\n', x0_val, y0_val);
fprintf('Antiderivative H(y) = F(x) : %s = %s + c\n', sol.separated_H_y, sol.separated_F_x);
fprintf('Calibrated Constant (c)    : %.6f\n', c_solved);
fprintf('Explicit Solution y(x)     : %s\n', sol.analytic_solution);
fprintf('Evaluation at x = %.4f     : Analytic = %.8f | ode45 = %.8f (Diff = %.2e)\n', ...
    x_target, y_analytic_target, y_num_target, discrepancy);
end