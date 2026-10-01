function report = solve_integral_completing_square(b_coeff, c_coeff, x_bounds)
% SOLVE_INTEGRAL_COMPLETING_SQUARE Solves definite integrals of the form:
% I = \int_{x1}^{x2} \frac{dx}{\sqrt{x^2 + b*x + c}} via completing the square (Lecture 28).
%
% Algebraic Formulation:
%   x^2 + b*x + c = (x + b/2)^2 + (c - b^2/4)
% Case (c - b^2/4 < 0, let d^2 = b^2/4 - c):
%   \int \frac{dv}{\sqrt{v^2 - d^2}} = ln|v + \sqrt{v^2 - d^2}| + C, where v = x + b/2.
%
% Syntax:
%   report = solve_integral_completing_square(b_coeff, c_coeff, x_bounds)
%
% Inputs:
%   b_coeff  - Linear coefficient b (default = 4.0 for x^2 + 4x)
%   c_coeff  - Constant term c (default = 0.0 for x^2 + 4x)
%   x_bounds - Two-element limits vector [x1, x2] in real valid domain (default: [1.0, 3.0])
%
% Outputs:
%   report   - Struct containing square completion shift, symbolic antiderivative, and values

arguments
    b_coeff  (1,1) double = 4.0
    c_coeff  (1,1) double = 0.0
    x_bounds (1,2) double = [1.0, 3.0]
end

x1 = x_bounds(1);
x2 = x_bounds(2);

if x1 >= x2
    error('Lower limit x1 must be strictly less than upper limit x2.');
end

% Verify real integrand domain (x^2 + b*x + c > 0 on [x1, x2])
test_poly = @(x) x.^2 + b_coeff * x + c_coeff;
if test_poly(x1) <= 0 || test_poly(x2) <= 0
    error('Domain violation: Argument under radical must be strictly positive on [x1, x2].');
end

% 1. Completing the Square: (x + shift)^2 + rem_const
shift_val = b_coeff / 2.0;
rem_const = c_coeff - (b_coeff^2) / 4.0;

% 2. Analytic Antiderivative Construction:
% For any x^2 + bx + c > 0, the universal antiderivative is:
% F(x) = ln | (x + b/2) + sqrt(x^2 + bx + c) |
F_anti = @(x) log(abs((x + shift_val) + sqrt(x.^2 + b_coeff * x + c_coeff)));

val_analytic = F_anti(x2) - F_anti(x1);

% 3. Numerical Adaptive Quadrature Verification
f_integrand = @(x) 1.0 ./ sqrt(x.^2 + b_coeff * x + c_coeff);
val_numerical = integral(f_integrand, x1, x2, 'RelTol', 1e-12, 'AbsTol', 1e-14);

discrepancy = abs(val_analytic - val_numerical);

% Assemble structured report
report.b_coefficient = b_coeff;
report.c_coefficient = c_coeff;
report.completed_square_form = sprintf('(x + %.4f)^2 + (%.4f)', shift_val, rem_const);
report.integration_limits = [x1, x2];
report.analytic_antiderivative = sprintf('ln|(x + %.4f) + \\surd(x^2 + %.2fx + %.2f)|', ...
    shift_val, b_coeff, c_coeff);
report.analytic_integral_value = val_analytic;
report.numerical_integral_value = val_numerical;
report.discrepancy = discrepancy;

% Console report
fprintf('=== Completing the Square Integral Solver ===\n');
fprintf('Integrand Radicand Form     : x^2 + %.2f*x + %.2f\n', b_coeff, c_coeff);
fprintf('Completed Square Formulation: %s\n', report.completed_square_form);
fprintf('Symbolic Antiderivative F(x): %s + C\n', report.analytic_antiderivative);
fprintf('Integration Interval [x1, x2]: [%.4f, %.4f]\n', x1, x2);
fprintf('Analytic Integral Value     : %.10f\n', val_analytic);
fprintf('Numerical Quadrature Check  : %.10f (Diff: %.2e)\n', val_numerical, discrepancy);
end