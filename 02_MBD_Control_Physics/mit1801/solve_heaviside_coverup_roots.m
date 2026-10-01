function sol = solve_heaviside_coverup_roots(numerator_poly, denominator_roots)
% SOLVE_HEAVISIDE_COVERUP_ROOTS Solves partial fraction decomposition for rational
% functions with distinct real roots: P(x) / \prod (x - p_i) = \Sum [ A_i / (x - p_i) ]
% Uses Heaviside's Cover-Up Method: A_k = P(p_k) / \prod_{j ~= k} (p_k - p_j).
%
% Syntax:
%   sol = solve_heaviside_coverup_roots(numerator_poly, denominator_roots)
%
% Inputs:
%   numerator_poly    - Polynomial coefficient vector P(x) in descending powers (default: [4, -1])
%   denominator_roots - Vector of distinct real roots [p_1, p_2, ...] (default: [1.0, -2.0])
%
% Outputs:
%   sol               - Struct containing computed residues, poles, and reconstructed verification

arguments
    numerator_poly    (1,:) double = [4.0, -1.0]
    denominator_roots (1,:) double = [1.0, -2.0]
end

num_poles = length(denominator_roots);
deg_num = length(numerator_poly) - 1;

if deg_num >= num_poles
    error('Degree of numerator P(x) must be strictly less than denominator Q(x) for standard PFE.');
end

% Verify roots are distinct
if length(unique(denominator_roots)) ~= num_poles
    error('Denominator roots must be strictly distinct for simple Heaviside cover-up.');
end

residues_A = zeros(num_poles, 1);

% 1. Heaviside Cover-Up Method Execution
for k = 1:num_poles
    pk = denominator_roots(k);
    % Evaluate Numerator at pk: P(pk)
    num_val = polyval(numerator_poly, pk);

    % Evaluate Denominator product excluding factor (x - pk)
    den_product = 1.0;
    for j = 1:num_poles
        if j ~= k
            den_product = den_product * (pk - denominator_roots(j));
        end
    end

    residues_A(k) = num_val / den_product;
end

% 2. Cross-Verification via Matrix System (Undetermined Coefficients)
% Reconstruct denominator polynomial Q(x) from roots
den_poly = poly(denominator_roots);
[matlab_r, matlab_p, ~] = residue(numerator_poly, den_poly);

% Sort for comparison alignment
[poles_sorted, idx_sort] = sort(denominator_roots');
res_sorted = residues_A(idx_sort);
[mat_p_sorted, idx_m] = sort(matlab_p);
mat_r_sorted = matlab_r(idx_m);

discrepancy = max(abs(res_sorted - mat_r_sorted));

% Construct symbolic antiderivative string
sym_anti_str = '';
for k = 1:num_poles
    sym_anti_str = sprintf('%s%+.4f*ln|x - (%.2f)| ', sym_anti_str, residues_A(k), denominator_roots(k));
end

% Assemble solution structure
sol.numerator = numerator_poly;
sol.poles = denominator_roots';
sol.residues_heaviside = residues_A;
sol.analytic_integral_string = [sym_anti_str, '+ c'];
sol.discrepancy_vs_residue = discrepancy;
sol.is_coverup_verified = (discrepancy < 1e-12);

% Console report
fprintf('=== Heaviside Cover-Up Partial Fraction Solver ===\n');
num_str = strjoin(arrayfun(@(c) sprintf('%+.4g', c), numerator_poly, 'UniformOutput', false), ' ');
fprintf('Numerator Polynomial P(x)  : %s\n', num_str);
fprintf('Denominator Factorized Q(x): \\prod (x - p_i) with %d Distinct Real Poles\n', num_poles);
fprintf('Computed Partial Fraction Decomposition:\n');
for k = 1:num_poles
    fprintf('  Pole %d at x = %+.4f | Residue A_%d = %+.6f\n', ...
        k, denominator_roots(k), k, residues_A(k));
end
fprintf('Symbolic Antiderivative F(x): %s\n', sol.analytic_integral_string);
fprintf('Residue Numerical Match Check: %s (Max Error: %.2e)\n', ...
    mat2str(sol.is_coverup_verified), discrepancy);
end