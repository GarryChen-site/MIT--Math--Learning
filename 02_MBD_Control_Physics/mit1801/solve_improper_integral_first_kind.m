function sol = solve_improper_integral_first_kind(f_sym, lower_bound_a, target_tolerance_eps)
% SOLVE_IMPROPER_INTEGRAL_FIRST_KIND Evaluates improper integrals of the first kind:
% I = \int_a^\infty f(x) dx = lim_{M -> \infty} \int_a^M f(x) dx (Lecture 35).
% Determines convergence, evaluates exact analytical limits, and computes the
% minimal simulation truncation horizon M_eps such that | \int_{M_eps}^\infty f(x) dx | < eps.
%
% Syntax:
%   sol = solve_improper_integral_first_kind(f_sym, lower_bound_a, target_tolerance_eps)
%
% Inputs:
%   f_sym                - Symbolic integrand expression in variable 'x' (default: exp(-0.5*x))
%   lower_bound_a        - Real finite lower boundary a (default = 0.0)
%   target_tolerance_eps - Precision truncation threshold for simulation horizon (default = 1e-4)
%
% Outputs:
%   sol                  - Struct containing convergence flag, exact value, M_eps, and numerical check

arguments
    f_sym                = []
    lower_bound_a        (1,1) double = 0.0
    target_tolerance_eps (1,1) double {mustBePositive} = 1e-4
end

syms x M real positive
if isempty(f_sym), f_sym = exp(-0.5 * x); end

% 1. Formulate Parametric Finite Definite Integral: F(M) = \int_a^M f(x) dx
F_M_sym = int(f_sym, x, lower_bound_a, M);

% 2. Evaluate Limit as M -> \infty
I_limit_sym = limit(F_M_sym, M, Inf);

is_convergent = ~isinf(double(I_limit_sym)) && ~isnan(double(I_limit_sym));

if is_convergent
    val_analytic = double(I_limit_sym);

    % 3. Truncation Horizon Evaluation: Residual tail \int_{M_eps}^\infty f(x) dx <= eps
    % Equivalent to solving: I_limit_sym - F_M_sym == target_tolerance_eps
    tail_sym = simplify(I_limit_sym - F_M_sym);
    try
        M_solved = double(solve(tail_sym == target_tolerance_eps, M, 'Real', true));
        valid_M = M_solved(M_solved > lower_bound_a);
        if ~isempty(valid_M)
            M_horizon = valid_M(1);
        else
            M_horizon = NaN;
        end
    catch
        M_horizon = NaN;
    end

    % Numerical Verification via Adaptive Quadrature with Inf limit
    f_handle = matlabFunction(f_sym, 'Vars', x);
    val_numerical = integral(f_handle, lower_bound_a, Inf, 'RelTol', 1e-12, 'AbsTol', 1e-14);
    discrepancy = abs(val_analytic - val_numerical);

else
    val_analytic = Inf;
    val_numerical = Inf;
    M_horizon = NaN;
    discrepancy = NaN;
end

% Assemble solution structure
sol.integrand = char(f_sym);
sol.lower_limit_a = lower_bound_a;
sol.tolerance_target = target_tolerance_eps;
sol.is_convergent = is_convergent;
sol.analytic_value = val_analytic;
sol.numerical_value = val_numerical;
sol.discrepancy = discrepancy;
sol.required_truncation_horizon_M = M_horizon;

% Console report
fprintf('=== Improper Integral (First Kind) Evaluator ===\n');
fprintf('Integrand f(x)              : %s\n', sol.integrand);
fprintf('Integration Interval        : [%.2f, +\\infty)\n', lower_bound_a);
fprintf('Convergence Status          : %s\n', mat2str(is_convergent));
if is_convergent
    fprintf('Exact Analytical Value      : %.10f\n', val_analytic);
    fprintf('Adaptive Quadrature Check   : %.10f (Diff: %.2e)\n', val_numerical, discrepancy);
    if ~isnan(M_horizon)
        fprintf('Truncation Horizon (M_eps)  : M >= %.4f s (Ensures Tail Residual < %.1e)\n', ...
            M_horizon, target_tolerance_eps);
    end
else
    fprintf('Integral Behavior           : Divergent to +Infinity\n');
end
end