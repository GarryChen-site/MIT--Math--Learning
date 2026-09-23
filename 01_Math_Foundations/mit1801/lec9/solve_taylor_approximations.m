function sol = solve_taylor_approximations(f_sym, x0_val, delta_x_val)
% SOLVE_TAYLOR_APPROXIMATIONS Computes 1st-order linear and 2nd-order quadratic
% Taylor approximations for an arbitrary function f(x) at expansion basepoint x0.
% Formulas:
%   Linear:    L(x) = f(x0) + f'(x0) * (x - x0)
%   Quadratic: Q(x) = f(x0) + f'(x0) * (x - x0) + 0.5 * f''(x0) * (x - x0)^2
%
% Syntax:
%   sol = solve_taylor_approximations(f_sym, x0_val, delta_x_val)
%
% Inputs:
%   f_sym       - Symbolic expression in variable 'x' (default: exp(-2*x)/sqrt(1+x))
%   x0_val      - Basepoint operating coordinate (default = 0.0)
%   delta_x_val - Perturbation displacement (x - x0) for accuracy check (default = 0.05)
%
% Outputs:
%   sol         - Struct containing coefficients, symbolic approximations, and errors

arguments
    f_sym       = []
    x0_val      (1,1) double = 0.0
    delta_x_val (1,1) double = 0.05
end

syms x
if isempty(f_sym), f_sym = exp(-2*x) / sqrt(1 + x); end

% 1. Successive symbolic derivatives
df_dx = diff(f_sym, x);
d2f_dx2 = diff(df_dx, x);

% 2. Evaluate derivatives at basepoint x0
f0 = double(subs(f_sym, x, x0_val));
df0 = double(subs(df_dx, x, x0_val));
d2f0 = double(subs(d2f_dx2, x, x0_val));

% 3. Formulate polynomial approximation functions
linear_poly_sym = simplify(f0 + df0 * (x - x0_val));
quadratic_poly_sym = simplify(linear_poly_sym + 0.5 * d2f0 * (x - x0_val)^2);

% 4. Numerical error validation at x = x0 + delta_x
x_test = x0_val + delta_x_val;
val_exact = double(subs(f_sym, x, x_test));
val_linear = double(subs(linear_poly_sym, x, x_test));
val_quadratic = double(subs(quadratic_poly_sym, x, x_test));

err_linear = abs(val_exact - val_linear);
err_quadratic = abs(val_exact - val_quadratic);

% Populate solution structure
sol.expression = char(f_sym);
sol.basepoint = x0_val;
sol.delta_x = delta_x_val;
sol.f0 = f0;
sol.df0 = df0;
sol.d2f0 = d2f0;
sol.linear_polynomial = char(linear_poly_sym);
sol.quadratic_polynomial = char(quadratic_poly_sym);
sol.exact_value = val_exact;
sol.linear_approx_value = val_linear;
sol.quadratic_approx_value = val_quadratic;
sol.linear_absolute_error = err_linear;
sol.quadratic_absolute_error = err_quadratic;

% Console output
fprintf('=== Taylor Polynomial Approximations Analysis ===\n');
fprintf('Target Function f(x)       : %s\n', sol.expression);
fprintf('Expansion Basepoint (x0)   : %.4f\n', x0_val);
fprintf('Derivatives at Basepoint   : f(x0)=%.4f, f''(x0)=%.4f, f''''(x0)=%.4f\n', f0, df0, d2f0);
fprintf('Linear Polynomial L(x)     : %s\n', sol.linear_polynomial);
fprintf('Quadratic Polynomial Q(x)  : %s\n', sol.quadratic_polynomial);
fprintf('Evaluation at x = %.4f (Delta x = %.4f):\n', x_test, delta_x_val);
fprintf('  Exact Value              : %.8f\n', val_exact);
fprintf('  Linear Approx Value      : %.8f (Absolute Error: %.2e)\n', val_linear, err_linear);
fprintf('  Quadratic Approx Value   : %.8f (Absolute Error: %.2e)\n', val_quadratic, err_quadratic);
end