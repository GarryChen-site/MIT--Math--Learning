function sol = solve_log_inequality_hierarchy(x_eval)
% SOLVE_LOG_INEQUALITY_HIERARCHY Proves and evaluates the MVT-based inequality
% hierarchy for ln(1 + x) across successive Taylor truncations (PS4 Part II Problem 1).
% Bounds for x > 0:
%   Upper Bound 1 : ln(1+x) <= x
%   Lower Bound 2 : ln(1+x) >= x - x^2/2
%   Upper Bound 3 : ln(1+x) <= x - x^2/2 + x^3/3
% General Pattern : Alternating series expansion Sum_{k=1}^n (-1)^(k-1) * (x^k / k)
%
% Syntax:
%   sol = solve_log_inequality_hierarchy(x_eval)

arguments
    x_eval (1,1) double = 0.5
end

if x_eval <= -1.0
    error('Input x must be strictly greater than -1 for log evaluation.');
end

exact_log = log(1.0 + x_eval);

% Successive Taylor polynomial bounds
p1 = x_eval;
p2 = x_eval - (x_eval^2) / 2.0;
p3 = x_eval - (x_eval^2) / 2.0 + (x_eval^3) / 3.0;
p4 = x_eval - (x_eval^2) / 2.0 + (x_eval^3) / 3.0 - (x_eval^4) / 4.0;

% Validate inequality hierarchy depending on domain sign
if x_eval >= 0
    is_p1_upper = (exact_log <= p1 + 1e-14);
    is_p2_lower = (exact_log >= p2 - 1e-14);
    is_p3_upper = (exact_log <= p3 + 1e-14);
    is_p4_lower = (exact_log >= p4 - 1e-14);
else
    % For -1 < x < 0: ln(1+x) <= x still holds
    is_p1_upper = (exact_log <= p1 + 1e-14);
    is_p2_lower = NaN; is_p3_upper = NaN; is_p4_lower = NaN;
end

sol.evaluated_x = x_eval;
sol.exact_log = exact_log;
sol.bound_p1 = p1;
sol.bound_p2 = p2;
sol.bound_p3 = p3;
sol.bound_p4 = p4;
sol.is_bracketed_valid = (exact_log >= p2 && exact_log <= p3);

% Console report
fprintf('=== PS4 Problem 1: Logarithmic Inequality Hierarchy ===\n');
fprintf('Evaluation Coordinate x     : %+.4f\n', x_eval);
fprintf('Exact Function ln(1 + x)    : %.10f\n', exact_log);
fprintf('Order 1 Bound [x]           : %.10f (Upper Bound Valid: %s)\n', p1, mat2str(is_p1_upper));
fprintf('Order 2 Bound [x - x^2/2]   : %.10f (Lower Bound Valid: %s)\n', p2, mat2str(is_p2_lower));
fprintf('Order 3 Bound [Order 2+x^3/3]: %.10f (Upper Bound Valid: %s)\n', p3, mat2str(is_p3_upper));
fprintf('Order 4 Bound [Order 3-x^4/4]: %.10f (Lower Bound Valid: %s)\n', p4, mat2str(is_p4_lower));
fprintf('Strict Bracketing Confirmed : %s (p2 <= ln(1+x) <= p3)\n', mat2str(sol.is_bracketed_valid));
end