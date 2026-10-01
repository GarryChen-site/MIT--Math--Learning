function report = solve_partial_fraction_irreducible_quad(alpha_pole, beta_sq, x_bounds)
% SOLVE_PARTIAL_FRACTION_IRREDUCIBLE_QUAD Solves partial fraction expansion
% for rational functions with an irreducible quadratic term (Lecture 29 Example 3):
% Form: f(x) = 1 / [ (x^2 + beta_sq) * (x - alpha_pole) ]
% Expansion: A / (x - alpha) + (B*x + C) / (x^2 + beta_sq)
% Integral: A*ln|x - alpha| + (B/2)*ln(x^2 + beta_sq) + (C/beta)*atan(x/beta) + const
%
% Syntax:
%   report = solve_partial_fraction_irreducible_quad(alpha_pole, beta_sq, x_bounds)
%
% Inputs:
%   alpha_pole - Real root of linear factor (default = 1.0 for x - 1)
%   beta_sq    - Positive constant in irreducible quadratic (default = 1.0 for x^2 + 1)
%   x_bounds   - Real limits [x1, x2] avoiding x = alpha (default: [2.0, 5.0])
%
% Outputs:
%   report     - Struct containing coefficients A, B, C, symbolic integral, and numerical check

arguments
    alpha_pole (1,1) double = 1.0
    beta_sq    (1,1) double {mustBePositive} = 1.0
    x_bounds   (1,2) double = [2.0, 5.0]
end

x1 = x_bounds(1);
x2 = x_bounds(2);

if x1 >= x2
    error('Lower limit x1 must be strictly less than upper limit x2.');
end

if (x1 <= alpha_pole && x2 >= alpha_pole)
    error('Integration interval [%.2f, %.2f] contains pole singularity x = %.2f.', x1, x2, alpha_pole);
end

% 1. Algebraic Derivation from Lecture 29:
% 1 = A*(x^2 + beta_sq) + (B*x + C)*(x - alpha)
% Cover-up for A at x = alpha:
A_coeff = 1.0 / (alpha_pole^2 + beta_sq);

% Plug in x = 0: 1 = A*beta_sq - C*alpha => C = (A*beta_sq - 1) / alpha
if abs(alpha_pole) > 1e-12
    C_coeff = (A_coeff * beta_sq - 1.0) / alpha_pole;
else
    % If alpha = 0, direct substitution
    C_coeff = 0.0;
end

% Equating x^2 coefficient: 0 = A + B => B = -A
B_coeff = -A_coeff;

beta_root = sqrt(beta_sq);

% 2. Antiderivative Function Handle:
% F(x) = A*ln|x - alpha| + (B/2)*ln(x^2 + beta_sq) + (C / beta)*atan(x / beta)
F_anti = @(x) A_coeff * log(abs(x - alpha_pole)) + ...
    (B_coeff / 2.0) * log(x.^2 + beta_sq) + ...
    (C_coeff / beta_root) * atan(x ./ beta_root);

analytic_val = F_anti(x2) - F_anti(x1);

% 3. Numerical Quadrature Verification
f_rational = @(x) 1.0 ./ ((x.^2 + beta_sq) .* (x - alpha_pole));
numeric_val = integral(f_rational, x1, x2, 'RelTol', 1e-12, 'AbsTol', 1e-14);

discrepancy = abs(analytic_val - numeric_val);

% Assemble structured report
report.alpha_pole = alpha_pole;
report.beta_squared = beta_sq;
report.coefficient_A = A_coeff;
report.coefficient_B = B_coeff;
report.coefficient_C = C_coeff;
report.partial_fraction_str = sprintf('%.4f/(x - %.2f) + (%.4f*x %+.4f)/(x^2 + %.2f)', ...
    A_coeff, alpha_pole, B_coeff, C_coeff, beta_sq);
report.integration_limits = [x1, x2];
report.analytic_integral_value = analytic_val;
report.numerical_integral_value = numeric_val;
report.discrepancy = discrepancy;

% Console report
fprintf('=== Irreducible Quadratic Partial Fraction Solver ===\n');
fprintf('Rational Plant Model       : 1 / [ (x^2 + %.2f) * (x - %.2f) ]\n', beta_sq, alpha_pole);
fprintf('Decomposed Partial Fraction: %s\n', report.partial_fraction_str);
fprintf('Resolved Coefficients      : A = %+.6f | B = %+.6f | C = %+.6f\n', ...
    A_coeff, B_coeff, C_coeff);
fprintf('Integration Window [x1, x2]: [%.4f, %.4f]\n', x1, x2);
fprintf('Analytic Integral Result   : %.10f\n', analytic_val);
fprintf('Numerical Quadrature Check : %.10f (Diff: %.2e)\n', numeric_val, discrepancy);
end