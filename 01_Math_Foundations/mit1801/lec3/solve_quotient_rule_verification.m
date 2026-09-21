function verification = solve_quotient_rule_verification(u_sym, v_sym, eval_point)
% SOLVE_QUOTIENT_RULE_VERIFICATION Compares analytical quotient rule derivative
% with high-order numerical finite difference approximation.
%
% Syntax:
%   verification = solve_quotient_rule_verification(u_sym, v_sym, eval_point)
%
% Inputs:
%   u_sym      - Symbolic expression for numerator u(x)
%   v_sym      - Symbolic expression for denominator v(x)
%   eval_point - Numeric scalar for point of evaluation
%
% Outputs:
%   verification - Struct containing analytical derivative, quotient rule steps, and error

arguments
    u_sym      = []
    v_sym      = []
    eval_point (1,1) double = 2.0
end

syms x
if isempty(u_sym), u_sym = x^2 + 1; end
if isempty(v_sym), v_sym = cos(x) + 2; end

% Symbolic derivatives of individual terms
du_sym = diff(u_sym, x);
dv_sym = diff(v_sym, x);

% Quotient rule formula: (u'v - uv') / v^2
quotient_formula_sym = (du_sym * v_sym - u_sym * dv_sym) / (v_sym^2);


% Evaluate numerically
val_analytical = double(subs(quotient_formula_sym, x, eval_point));


% Central finite difference verification
h_step = 1e-6;
f_num = @(val) double(subs(u_sym / v_sym, x, val));
val_numerical = (f_num(eval_point + h_step) - f_num(eval_point - h_step)) / (2 * h_step);

abs_error = abs(val_analytical - val_numerical);

% Populate output structure
verification.eval_point = eval_point;
verification.quotient_expression = char(simplify(quotient_formula_sym));
verification.analytical_value = val_analytical;
verification.numerical_value = val_numerical;
verification.absolute_difference = abs_error;

% Console output
fprintf('=== Quotient Rule Derivative Verification ===\n');
fprintf('Rational Function          : (%s) / (%s)\n', char(u_sym), char(v_sym));
fprintf('Analytical Formula (u/v)''  : %s\n', verification.quotient_expression);
fprintf('Evaluated at x = %.4f      : Analytical = %.8f | Numerical = %.8f\n', ...
    eval_point, val_analytical, val_numerical);
fprintf('Absolute Discrepancy       : %.2e\n', abs_error);
end