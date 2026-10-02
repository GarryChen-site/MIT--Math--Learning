function report = evaluate_power_indeterminate_form(base_f_sym, exponent_g_sym, limit_point)
% EVALUATE_POWER_INDETERMINATE_FORM Evaluates indeterminate power limits:
% L = lim_{x -> a} [f(x)]^[g(x)] (forms: 0^0, 1^inf, inf^0) using logarithmic transformation:
% ln(L) = lim_{x -> a} [ ln(f(x)) / (1 / g(x)) ] via L'Hopital's rule => L = exp(ln(L)).
%
% Syntax:
%   report = evaluate_power_indeterminate_form(base_f_sym, exponent_g_sym, limit_point)
%
% Inputs:
%   base_f_sym     - Symbolic base expression f(x) (default: x for x^x)
%   exponent_g_sym - Symbolic exponent expression g(x) (default: x for x^x)
%   limit_point    - Limit coordinate a (default = 0.0)
%
% Outputs:
%   report         - Struct containing exponent limit, total limit, and classification

arguments
    base_f_sym     = []
    exponent_g_sym = []
    limit_point    (1,1) double = 0.0
end

syms x real
if isempty(base_f_sym), base_f_sym = x; end
if isempty(exponent_g_sym), exponent_g_sym = x; end

% 1. Classify indeterminate power structure at target point
f_at_a = limit(base_f_sym, x, limit_point, 'Right');
g_at_a = limit(exponent_g_sym, x, limit_point, 'Right');

f_val = double(f_at_a);
g_val = double(g_at_a);

if abs(f_val) < 1e-12 && abs(g_val) < 1e-12
    indeterminate_type = 'Type 0^0 (Zero to the Zero Power)';
elseif abs(f_val - 1.0) < 1e-12 && isinf(g_val)
    indeterminate_type = 'Type 1^\infty (Unit to Infinite Power)';
elseif isinf(f_val) && abs(g_val) < 1e-12
    indeterminate_type = 'Type \infty^0 (Infinity to Zero Power)';
else
    indeterminate_type = 'Determinate or Standard Form';
end

% 2. Transform to Fraction for L'Hopital: ln(L) = lim [ ln(f(x)) / (1 / g(x)) ]
ln_numerator = log(base_f_sym);
inv_denominator = 1.0 / exponent_g_sym;

% Solve transformed limit using Symbolic Math Engine
ln_L_sym = limit(ln_numerator / inv_denominator, x, limit_point, 'Right');
ln_L_val = double(ln_L_sym);

% Exponentiate to get final limit L = exp(ln(L))
L_final = exp(ln_L_val);

% 3. Numerical Verification via Small Perturbation from the Right
h = 1e-6;
x_test = limit_point + h;
f_num = matlabFunction(base_f_sym, 'Vars', x);
g_num = matlabFunction(exponent_g_sym, 'Vars', x);
numerical_val = (f_num(x_test))^(g_num(x_test));

discrepancy = abs(L_final - numerical_val);

% Assemble structured report
report.base_function = char(base_f_sym);
report.exponent_function = char(exponent_g_sym);
report.limit_target = limit_point;
report.indeterminate_classification = indeterminate_type;
report.log_transformed_limit = ln_L_val;
report.analytic_limit_value = L_final;
report.numerical_perturbation_value = numerical_val;
report.discrepancy = discrepancy;

% Console report
fprintf('=== Indeterminate Power Form Evaluator [f(x)^g(x)] ===\n');
fprintf('Target Limit Expression    : lim_{x -> %+.4f+} [ %s ]^[ %s ]\n', ...
    limit_point, report.base_function, report.exponent_function);
fprintf('Indeterminate Structure    : %s\n', indeterminate_type);
fprintf('Transformed Exponent Limit : ln(L) = %.10f\n', ln_L_val);
fprintf('Final Analytic Limit (L)   : %.10f (e^{ln(L)})\n', L_final);
fprintf('Numerical Check at x = a+h : %.10f (Diff: %.2e)\n', numerical_val, discrepancy);
end