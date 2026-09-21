function result = solve_higher_order_derivatives(expr_sym, max_order, eval_point)
% SOLVE_HIGHER_ORDER_DERIVATIVES Computes successive derivatives up to order n,
% validates factorial property D^n(x^n) = n!, and evaluates states numerically.
%
% Syntax:
%   result = solve_higher_order_derivatives(expr_sym, max_order, eval_point)
%
% Inputs:
%   expr_sym   - Symbolic expression in variable 'x'
%   max_order  - Maximum derivative order (integer >= 1, default = 4)
%   eval_point - Operating point x0 for numeric evaluation
%
% Outputs:
%   result     - Struct containing symbolic derivative table, values, and factorial check

arguments
    expr_sym   = []
    max_order  (1,1) double {mustBeInteger, mustBePositive} = 4
    eval_point (1,1) double = 1.0
end

syms x
if isempty(expr_sym), expr_sym = x^max_order; end

deriv_sym = cell(max_order, 1);
deriv_val = zeros(max_order, 1);

current_deriv = expr_sym;
for k = 1:max_order
    current_deriv = diff(current_deriv, x);
    deriv_sym{k} = simplify(current_deriv);
    deriv_val(k) = double(subs(current_deriv, x, eval_point));
end

% Factorial check for monomial x^n: D^n(x^n) == factorial(n)
expected_factorial = factorial(max_order);
is_monomial_case = isequal(expr_sym, x^max_order);
factorial_match = is_monomial_case && (deriv_val(max_order) == expected_factorial);

result.expression = char(expr_sym);
result.order = max_order;
result.eval_point = eval_point;
result.symbolic_derivatives = deriv_sym;
result.evaluated_values = deriv_val;
result.is_factorial_verified = factorial_match;

% Console output
fprintf('=== Higher-Order Derivatives Analysis (Up to Order %d) ===\n', max_order);
fprintf('Base Function f(x) : %s\n', result.expression);
fprintf('Operating Point x0 : %.4f\n', eval_point);
for k = 1:max_order
    fprintf('  Order %d [D^%d f]   : %s  ==> Value: %.4f\n', ...
        k, k, char(deriv_sym{k}), deriv_val(k));
end
if is_monomial_case
    fprintf('Factorial Property Check: D^%d(x^%d) = %d (Theory: %d! = %d) -> Verified: %s\n', ...
        max_order, max_order, round(deriv_val(max_order)), max_order, expected_factorial, ...
        mat2str(factorial_match));
end
end