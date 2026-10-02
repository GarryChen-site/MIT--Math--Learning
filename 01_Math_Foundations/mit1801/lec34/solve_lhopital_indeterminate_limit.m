function sol = solve_lhopital_indeterminate_limit(f_sym, g_sym, limit_point, max_depth)
% SOLVE_LHOPITAL_INDETERMINATE_LIMIT Evaluates lim_{x -> a} [f(x) / g(x)]
% using iterative L'Hopital's rule with strict pre-condition verification
% ("Look before you L'Hop" guard to avoid false zero/infinite convergence).
%
% Syntax:
%   sol = solve_lhopital_indeterminate_limit(f_sym, g_sym, limit_point, max_depth)
%
% Inputs:
%   f_sym       - Symbolic numerator in variable 'x' (default: cos(x) - 1)
%   g_sym       - Symbolic denominator in variable 'x' (default: x^2)
%   limit_point - Limit target coordinate a (scalar or +/-Inf, default = 0.0)
%   max_depth   - Maximum allowed recursive L'Hopital differentiations (default = 5)
%
% Outputs:
%   sol         - Struct containing limit value, number of differentiations, and status

arguments
    f_sym       = []
    g_sym       = []
    limit_point (1,1) double = 0.0
    max_depth   (1,1) double {mustBeInteger, mustBePositive} = 5
end

syms x real
if isempty(f_sym), f_sym = cos(x) - 1; end
if isempty(g_sym), g_sym = x^2; end

curr_f = f_sym;
curr_g = g_sym;

derivation_steps = cell(max_depth + 1, 3);
depth_count = 0;
resolved_limit = NaN;

for k = 1:max_depth
    % 1. Evaluate current limits of numerator and denominator
    f_val = limit(curr_f, x, limit_point);
    g_val = limit(curr_g, x, limit_point);

    derivation_steps{k, 1} = char(curr_f);
    derivation_steps{k, 2} = char(curr_g);
    derivation_steps{k, 3} = sprintf('f(a)=%s, g(a)=%s', char(f_val), char(g_val));

    % Check if indeterminate form 0/0 or Inf/Inf exists
    is_zero_zero = (abs(double(f_val)) < 1e-12) && (abs(double(g_val)) < 1e-12);
    is_inf_inf = (isinf(double(f_val))) && (isinf(double(g_val)));

    if is_zero_zero || is_inf_inf
        % Apply L'Hopital rule: differentiate both top and bottom
        depth_count = depth_count + 1;
        curr_f = diff(curr_f, x);
        curr_g = diff(curr_g, x);
    else
        % Form is determinate! Check denominator
        if abs(double(g_val)) > 1e-12
            resolved_limit = double(f_val / g_val);
        else
            % Denominator is 0 but numerator is non-zero: Essential singularity (Inf)
            resolved_limit = double(f_val / g_val);
        end
        break;
    end
end

% Final check after loop if terminated at max_depth
if isnan(resolved_limit)
    f_val = limit(curr_f, x, limit_point);
    g_val = limit(curr_g, x, limit_point);
    resolved_limit = double(f_val / g_val);
end

% Numerical perturbation validation
h = 1e-5;
if isinf(limit_point)
    x_num = 1e6;
else
    x_num = limit_point + h;
end
f_h = matlabFunction(f_sym, 'Vars', x);
g_h = matlabFunction(g_sym, 'Vars', x);
num_ratio_val = f_h(x_num) / g_h(x_num);

discrepancy = abs(resolved_limit - num_ratio_val);

% Assemble solution structure
sol.numerator_expr = char(f_sym);
sol.denominator_expr = char(g_sym);
sol.target_point = limit_point;
sol.lhopital_iterations_taken = depth_count;
sol.analytic_limit_value = resolved_limit;
sol.numerical_perturbation_value = num_ratio_val;
sol.discrepancy = discrepancy;
sol.is_determinate_resolved = ~isnan(resolved_limit);

% Console report
fprintf('=== L''Hopital Indeterminate Form Limit Solver ===\n');
fprintf('Rational Quotient          : [ %s ] / [ %s ]\n', sol.numerator_expr, sol.denominator_expr);
fprintf('Limit Target Point (a)     : x -> %.4f\n', limit_point);
fprintf('L''Hopital Operations Count : %d Differentials Applied\n', depth_count);
fprintf('Analytic Resolved Limit    : %.10f\n', resolved_limit);
fprintf('Numerical Perturbation Val : %.10f (Diff: %.2e)\n', num_ratio_val, discrepancy);
fprintf('Convergence Validity Check : %s\n', mat2str(sol.is_determinate_resolved));
end