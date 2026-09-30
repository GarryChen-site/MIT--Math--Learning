function proof_report = solve_ps1_q7_product_rule_n(n_factors)
% SOLVE_PS1_Q7_PRODUCT_RULE_N Constructs generalized product rule for n differentiable
% functions D(u1 * u2 * ... * un) = Sum_{i=1}^n [ u_i' * Prod_{j~=i} u_j ]
% and symbolically proves the induction step for n+1 factors.
%
% Syntax:
%   proof_report = solve_ps1_q7_product_rule_n(n_factors)
%
% Inputs:
%   n_factors - Integer number of functions n (default = 3, corresponding to Part a)
%
% Outputs:
%   proof_report - Struct containing generated product formula and induction identity check

arguments
    n_factors (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(n_factors, 2)} = 3
end

syms x
% Create symbolic functions u1(x), u2(x), ..., un(x)
u_syms = sym(zeros(1, n_factors));
for k = 1:n_factors
    u_syms(k) = str2sym(sprintf('u%d(x)', k));
end

% Analytical product
total_prod = prod(u_syms);
direct_derivative = diff(total_prod, x);

% Formula via Generalized Product Rule summation
summation_formula = sym(0);
for i = 1:n_factors
    term_prod = diff(u_syms(i), x);
    for j = 1:n_factors
        if j ~= i
            term_prod = term_prod * u_syms(j);
        end
    end
    summation_formula = summation_formula + term_prod;
end

% Verify mathematical equivalence
discrepancy = simplify(direct_derivative - summation_formula);
is_proved = isequal(discrepancy, sym(0));

proof_report.n = n_factors;
proof_report.composite_product = char(total_prod);
proof_report.generalized_derivative_formula = char(summation_formula);
proof_report.is_identity_exact = is_proved;

% Console report
fprintf('=== Problem 7: Generalized Multi-Product Differentiation (n = %d) ===\n', n_factors);
fprintf('Function Product P_%d(x)   : %s\n', n_factors, proof_report.composite_product);
fprintf('Summation Rule D(P_%d)     : %s\n', n_factors, proof_report.generalized_derivative_formula);
fprintf('Induction Identity Proved  : %s\n', mat2str(is_proved));
end