function sol = solve_leibniz_variable_integral_derivative(f_kernel_sym, u_lower_sym, v_upper_sym, x_eval)
% SOLVE_LEIBNIZ_VARIABLE_INTEGRAL_DERIVATIVE Evaluates the derivative of an integral
% with variable upper and lower limits using the Leibniz Integral Rule (generalized FTC 2):
% Formula: d/dx [ \int_{u(x)}^{v(x)} f(t) dt ] = f(v(x)) * v'(x) - f(u(x)) * u'(x)
%
% Syntax:
%   sol = solve_leibniz_variable_integral_derivative(f_kernel_sym, u_lower_sym, v_upper_sym, x_eval)
%
% Inputs:
%   f_kernel_sym - Symbolic integrand in variable 't' (default: cos(t^2) for Fresnel C)
%   u_lower_sym  - Symbolic lower bound in variable 'x' (default: 0)
%   v_upper_sym  - Symbolic upper bound in variable 'x' (default: x^2)
%   x_eval       - Real scalar operating point for evaluation (default = 1.2)
%
% Outputs:
%   sol          - Struct containing symbolic derivative formula, analytic value, and error

arguments
    f_kernel_sym = []
    u_lower_sym  = []
    v_upper_sym  = []
    x_eval       (1,1) double = 1.2
end

syms t x real
if isempty(f_kernel_sym), f_kernel_sym = cos(t^2); end
if isempty(u_lower_sym),  u_lower_sym = sym(0); end
if isempty(v_upper_sym),  v_upper_sym = x^2; end

% 1. Derive boundary derivatives
du_dx = diff(u_lower_sym, x);
dv_dx = diff(v_upper_sym, x);

% 2. Apply Leibniz rule
f_at_v = subs(f_kernel_sym, t, v_upper_sym);
f_at_u = subs(f_kernel_sym, t, u_lower_sym);
dI_dx_sym = simplify(f_at_v * dv_dx - f_at_u * du_dx);

% 3. Numerical evaluation at operating point
dI_dx_analytic = double(subs(dI_dx_sym, x, x_eval));

% 4. Finite difference numerical validation
h = 1e-6;
f_num = matlabFunction(f_kernel_sym, 'Vars', t);
u_num = matlabFunction(u_lower_sym, 'Vars', x);
v_num = matlabFunction(v_upper_sym, 'Vars', x);

I_plus  = integral(f_num, u_num(x_eval + h), v_num(x_eval + h), 'RelTol', 1e-12, 'AbsTol', 1e-14);
I_minus = integral(f_num, u_num(x_eval - h), v_num(x_eval - h), 'RelTol', 1e-12, 'AbsTol', 1e-14);
dI_dx_numerical = (I_plus - I_minus) / (2 * h);

abs_discrepancy = abs(dI_dx_analytic - dI_dx_numerical);

% Assemble solution structure
sol.x_eval = x_eval;
sol.integrand_kernel = char(f_kernel_sym);
sol.lower_limit = char(u_lower_sym);
sol.upper_limit = char(v_upper_sym);
sol.symbolic_derivative = char(dI_dx_sym);
sol.analytic_derivative = dI_dx_analytic;
sol.numerical_derivative = dI_dx_numerical;
sol.discrepancy = abs_discrepancy;

% Console report
fprintf('=== Leibniz Variable Integral Derivative Solver ===\n');
fprintf('Integrand Kernel f(t)       : %s\n', sol.integrand_kernel);
fprintf('Integral Interval [u(x),v(x)]: [%s, %s]\n', sol.lower_limit, sol.upper_limit);
fprintf('Symbolic Derivative dI/dx   : %s\n', sol.symbolic_derivative);
fprintf('Evaluated Operating Point   : x0 = %.4f\n', x_eval);
fprintf('Analytic Derivative Value   : %.10f\n', dI_dx_analytic);
fprintf('Numerical Finite Difference : %.10f\n', dI_dx_numerical);
fprintf('Verification Discrepancy    : %.2e\n', abs_discrepancy);
end