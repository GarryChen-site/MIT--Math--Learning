function sol = solve_trig_integral_reduction(n_pow, m_pow, x_bounds)
% SOLVE_TRIG_INTEGRAL_REDUCTION Solves indefinite and definite integrals of the form:
% \int sin^n(x) * cos^m(x) dx using Method A (odd substitution) or Method B (double-angle reduction).
%
% Syntax:
%   sol = solve_trig_integral_reduction(n_pow, m_pow, x_bounds)
%
% Inputs:
%   n_pow    - Integer power of sine (default = 3)
%   m_pow    - Integer power of cosine (default = 2)
%   x_bounds - Two-element limits vector [x1, x2] (default: [0.0, pi/2])
%
% Outputs:
%   sol      - Struct containing method classification, symbolic antiderivative, and values

arguments
    n_pow    (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    m_pow    (1,1) double {mustBeInteger, mustBeNonnegative} = 2
    x_bounds (1,2) double = [0.0, pi/2]
end

syms x real
integrand_sym = (sin(x)^n_pow) * (cos(x)^m_pow);

% 1. Strategy Classification from Lecture 26
if mod(n_pow, 2) == 1 || mod(m_pow, 2) == 1
    method_type = 'Method A: Odd power substitution (u = cos(x) or u = sin(x))';
else
    method_type = 'Method B: Even powers (Repeated double-angle reduction)';
end

% 2. Symbolic Antiderivative Integration
F_sym = int(integrand_sym, x);

% 3. Definite Integral Evaluation
x1 = x_bounds(1);
x2 = x_bounds(2);
val_analytic = double(subs(F_sym, x, x2) - subs(F_sym, x, x1));

% 4. High-Precision Numerical Validation
integrand_handle = @(xv) (sin(xv).^n_pow) .* (cos(xv).^m_pow);
val_numerical = integral(integrand_handle, x1, x2, 'RelTol', 1e-12, 'AbsTol', 1e-14);

discrepancy = abs(val_analytic - val_numerical);

% Assemble solution structure
sol.n_power = n_pow;
sol.m_power = m_pow;
sol.method_classification = method_type;
sol.integrand_expression = char(integrand_sym);
sol.antiderivative_expression = char(simplify(F_sym));
sol.bounds = [x1, x2];
sol.analytic_value = val_analytic;
sol.numerical_value = val_numerical;
sol.discrepancy = discrepancy;

% Console report
fprintf('=== Trigonometric Integral Reduction Solver ===\n');
fprintf('Powers                     : sin^%d(x) * cos^%d(x)\n', n_pow, m_pow);
fprintf('Recommended Strategy       : %s\n', method_type);
fprintf('Symbolic Antiderivative    : %s + c\n', sol.antiderivative_expression);
fprintf('Definite Limits [x1, x2]   : [%.4f, %.4f]\n', x1, x2);
fprintf('Analytic Integral Value    : %.10f\n', val_analytic);
fprintf('Numerical Quadrature Check : %.10f (Diff: %.2e)\n', val_numerical, discrepancy);
end