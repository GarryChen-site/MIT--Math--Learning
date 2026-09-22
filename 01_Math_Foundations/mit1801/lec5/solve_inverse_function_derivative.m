function result = solve_inverse_function_derivative(f_sym, x0_val)
% SOLVE_INVERSE_FUNCTION_DERIVATIVE Evaluates derivative of inverse function g(y) = f^(-1)(y)
% using the reciprocal derivative relationship: g'(y0) = 1 / f'(x0), where y0 = f(x0).
%
% Syntax:
%   result = solve_inverse_function_derivative(f_sym, x0_val)
%
% Inputs:
%   f_sym  - Symbolic forward function f(x)
%   x0_val - Evaluation point x0 for forward function
%
% Outputs:
%   result - Struct containing forward and inverse points, derivatives, and reciprocal check

arguments
    f_sym  = []
    x0_val (1,1) double = pi/4
end

syms x
if isempty(f_sym), f_sym = tan(x); end

% 1. Forward function and derivative
df_dx_sym = diff(f_sym, x);
y0_val = double(subs(f_sym, x, x0_val));
df_dx_val = double(subs(df_dx_sym, x, x0_val));

if abs(df_dx_val) < 1e-12
    error('Forward derivative is zero: Inverse derivative does not exist (Singularity).');
end

% 2. Inverse function derivative via reciprocal formula
dg_dy_reciprocal = 1.0 / df_dx_val;

% 3. Finite-difference numerical validation on inverse mapping
% We perturb y by dy and compute corresponding dx numerically
dy_step = 1e-6;
f_num = matlabFunction(f_sym);
% Numerical inversion via fzero
x_pert_plus  = fzero(@(xv) f_num(xv) - (y0_val + dy_step), x0_val);
x_pert_minus = fzero(@(xv) f_num(xv) - (y0_val - dy_step), x0_val);
dg_dy_numerical = (x_pert_plus - x_pert_minus) / (2 * dy_step);

abs_error = abs(dg_dy_reciprocal - dg_dy_numerical);

% Assemble results
result.f_symbolic = char(f_sym);
result.x0 = x0_val;
result.y0 = y0_val;
result.df_dx_at_x0 = df_dx_val;
result.dg_dy_at_y0 = dg_dy_reciprocal;
result.dg_dy_numerical = dg_dy_numerical;
result.discrepancy = abs_error;

% Console output
fprintf('=== Inverse Function Derivative Calculation ===\n');
fprintf('Forward Function f(x)         : %s\n', result.f_symbolic);
fprintf('Forward Operating Point (x0)  : %.4f\n', x0_val);
fprintf('Mapped Inverse Point (y0)     : %.4f\n', y0_val);
fprintf('Forward Derivative f''(x0)     : %.6f\n', df_dx_val);
fprintf('Inverse Derivative g''(y0)     : %.6f (1 / f''(x0))\n', dg_dy_reciprocal);
fprintf('Numerical Validation g''(y0)   : %.6f (Error: %.2e)\n', dg_dy_numerical, abs_error);
end