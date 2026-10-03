function sol = solve_taylor_polynomial_approximation(f_sym, base_b, target_order_N, test_interval)
% SOLVE_TAYLOR_POLYNOMIAL_APPROXIMATION Generates symbolic and numerical Taylor series
% polynomials for function f(x) expanded around operating base point b (Lecture 37).
% Formula: P_N(x) = \sum_{n=0}^N [ f^(n)(b) / n! ] * (x - b)^n
%
% Syntax:
%   sol = solve_taylor_polynomial_approximation(f_sym, base_b, target_order_N, test_interval)
%
% Inputs:
%   f_sym          - Symbolic function in variable 'x' (default: cos(x))
%   base_b         - Expansion operating point b (default = 0.0 for Maclaurin)
%   target_order_N - Maximum polynomial truncation order N (default = 4)
%   test_interval  - Two-element testing domain [x1, x2] for residual audit (default: [-pi, pi])
%
% Outputs:
%   sol            - Struct containing polynomial string, coefficients, and max truncation error

arguments
    f_sym          = []
    base_b         (1,1) double = 0.0
    target_order_N (1,1) double {mustBeInteger, mustBeNonnegative} = 4
    test_interval  (1,2) double = [-pi, pi]
end

syms x real
if isempty(f_sym), f_sym = cos(x); end

x1 = test_interval(1);
x2 = test_interval(2);

if x1 >= x2
    error('Invalid interval: Lower test bound must be strictly less than upper bound.');
end

% 1. Successive Derivative Calculation & Taylor Coefficients a_n = f^(n)(b) / n!
a_coeffs = zeros(target_order_N + 1, 1);
P_N_sym = sym(0);
current_deriv = f_sym;

for n = 0:target_order_N
    if n == 0
        f_n_at_b = double(subs(f_sym, x, base_b));
    else
        current_deriv = diff(current_deriv, x);
        f_n_at_b = double(subs(current_deriv, x, base_b));
    end

    coeff_val = f_n_at_b / factorial(n);
    a_coeffs(n+1) = coeff_val;

    % Build symbolic polynomial P_N(x)
    P_N_sym = P_N_sym + coeff_val * (x - base_b)^n;
end

% 2. Maximum Truncation Residual Evaluation on test_interval
f_handle = matlabFunction(f_sym, 'Vars', x);
P_handle = matlabFunction(P_N_sym, 'Vars', x);

x_test_grid = linspace(x1, x2, 1000);
y_true_grid = f_handle(x_test_grid);
y_poly_grid = P_handle(x_test_grid);

residuals = abs(y_true_grid - y_poly_grid);
[max_err, max_idx] = max(residuals);
worst_point = x_test_grid(max_idx);

% Assemble solution structure
sol.function_expression = char(f_sym);
sol.expansion_basepoint = base_b;
sol.truncation_order = target_order_N;
sol.taylor_coefficients = a_coeffs;
sol.symbolic_polynomial = char(simplify(P_N_sym));
sol.testing_interval = [x1, x2];
sol.max_truncation_error = max_err;
sol.worst_case_x = worst_point;

% Console report
fprintf('=== Taylor Series Polynomial Expansion Solver ===\n');
fprintf('Target Function f(x)       : %s\n', sol.function_expression);
fprintf('Operating Base Point (b)   : %.4f (Maclaurin if b = 0)\n', base_b);
fprintf('Truncation Polynomial Order: N = %d\n', target_order_N);
fprintf('Symbolic Approximation P(x): %s\n', sol.symbolic_polynomial);
fprintf('Taylor Coefficients List   : %s\n', mat2str(a_coeffs', 5));
fprintf('Max Truncation Error       : %.2e (Occurs at x = %+.4f on [%.2f, %.2f])\n', ...
    max_err, worst_point, x1, x2);
end