function result = solve_definite_integral_substitution(f_outer_sym, g_inner_sym, x_bounds)
% SOLVE_DEFINITE_INTEGRAL_SUBSTITUTION Solves definite integrals using u-substitution:
% Integral: \int_x1^x2 f(g(x)) * g'(x) dx = \int_{u1}^{u2} f(u) du
% Automatically remaps limits: u1 = g(x1), u2 = g(x2).
%
% Syntax:
%   result = solve_definite_integral_substitution(f_outer_sym, g_inner_sym, x_bounds)
%
% Inputs:
%   f_outer_sym - Outer integrand symbolic expression in variable 'u' (default: u^4)
%   g_inner_sym - Inner substitution expression in variable 'x' (default: x^3 + 2)
%   x_bounds    - Two-element limits vector [x1, x2] (default: [1.0, 2.0])
%
% Outputs:
%   result      - Struct containing remapped limits, analytic values, and numerical checks

arguments
    f_outer_sym = []
    g_inner_sym = []
    x_bounds    (1,2) double = [1.0, 2.0]
end

syms x u real
if isempty(f_outer_sym), f_outer_sym = u^4; end
if isempty(g_inner_sym), g_inner_sym = x^3 + 2; end

x1 = x_bounds(1);
x2 = x_bounds(2);

% 1. Compute inner derivative g'(x)
dg_dx = diff(g_inner_sym, x);

% 2. Full integrand in original domain x: f(g(x)) * g'(x)
integrand_x_sym = subs(f_outer_sym, u, g_inner_sym) * dg_dx;

% 3. Map integration limits to u-domain
u1 = double(subs(g_inner_sym, x, x1));
u2 = double(subs(g_inner_sym, x, x2));

% 4. Analytic evaluation in u-domain
F_u_sym = int(f_outer_sym, u);
val_analytic = double(subs(F_u_sym, u, u2) - subs(F_u_sym, u, u1));

% 5. High-precision numerical validation across both domains
handle_x = matlabFunction(integrand_x_sym);
handle_u = matlabFunction(f_outer_sym);

num_val_x = integral(handle_x, x1, x2, 'RelTol', 1e-12, 'AbsTol', 1e-14);
num_val_u = integral(handle_u, u1, u2, 'RelTol', 1e-12, 'AbsTol', 1e-14);

discrepancy = abs(val_analytic - num_val_x);

% Assemble result structure
result.original_bounds_x = [x1, x2];
result.substituted_bounds_u = [u1, u2];
result.inner_substitution_u = char(g_inner_sym);
result.outer_function_f = char(f_outer_sym);
result.original_integrand_x = char(simplify(integrand_x_sym));
result.analytic_integral_value = val_analytic;
result.numerical_integral_x = num_val_x;
result.numerical_integral_u = num_val_u;
result.validation_discrepancy = discrepancy;

% Console report
fprintf('=== Definite Integral u-Substitution Solver ===\n');
fprintf('Inner Substitution u(x)     : %s\n', result.inner_substitution_u);
fprintf('Transformed Integrand f(u)  : %s\n', result.outer_function_f);
fprintf('Original Bounds [x1, x2]    : [%.4f, %.4f]\n', x1, x2);
fprintf('Mapped Bounds [u1, u2]      : [%.4f, %.4f]\n', u1, u2);
fprintf('Analytic Integral Value     : %.10f\n', val_analytic);
fprintf('Original Domain Quadrature  : %.10f\n', num_val_x);
fprintf('Transformed Domain Quadrature: %.10f\n', num_val_u);
fprintf('Numerical Discrepancy Error : %.2e\n', discrepancy);
end