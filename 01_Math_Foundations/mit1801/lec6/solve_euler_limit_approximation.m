function result = solve_euler_limit_approximation(target_tol_pct)
% SOLVE_EULER_LIMIT_APPROXIMATION Solves Euler limit convergence (1 + 1/k)^k -> e
% and computes the minimum discretization subdivisions k required to meet a relative tolerance.
%
% Syntax:
%   result = solve_euler_limit_approximation(target_tol_pct)
%
% Inputs:
%   target_tol_pct - Target relative error percentage (e.g., 0.1 for 0.1% tolerance)
%
% Outputs:
%   result         - Struct containing required subdivisions, approximated e, and true error

arguments
    target_tol_pct (1,1) double {mustBePositive} = 0.5 % Default 0.5% error limit
end

e_exact = exp(1);
tol_ratio = target_tol_pct / 100.0;

% Asymptotic series expansion: (1 + 1/k)^k = e * (1 - 1/(2*k) + O(1/k^2))
% Relative error approx: |e - e_approx| / e approx 1 / (2*k) <= tol_ratio
% Hence: k_est >= 1 / (2 * tol_ratio)
k_initial_est = ceil(1.0 / (2.0 * tol_ratio));

% Precise numerical search using integer refinement
k_curr = max(1, k_initial_est - 10);
while true
    current_val = (1.0 + 1.0 / k_curr)^k_curr;
    current_err = abs(e_exact - current_val) / e_exact;
    if current_err <= tol_ratio
        k_min = k_curr;
        e_approx_val = current_val;
        final_err_pct = current_err * 100.0;
        break;
    end
    k_curr = k_curr + 1;
end

% Assemble result structure
result.target_tolerance_pct = target_tol_pct;
result.min_subdivisions_k = k_min;
result.approximated_e = e_approx_val;
result.exact_e = e_exact;
result.achieved_error_pct = final_err_pct;

% Console report
fprintf('=== Euler Limit Discretization Analysis ===\n');
fprintf('Target Relative Tolerance : %.4f%%\n', target_tol_pct);
fprintf('Required Minimum Steps (k): %d\n', k_min);
fprintf('Approximated Constant (e) : %.8f\n', e_approx_val);
fprintf('True Mathematical Base (e): %.8f\n', e_exact);
fprintf('Achieved Relative Error   : %.6f%%\n', final_err_pct);
end