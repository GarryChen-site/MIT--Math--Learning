function report = evaluate_p_integral_convergence(p_power, integral_kind, boundary_val)
% EVALUATE_P_INTEGRAL_CONVERGENCE Analyzes canonical p-power improper integrals:
%   Type 1 (Infinite Interval): \int_a^\infty (1 / x^p) dx  (Converges if p > 1, Diverges if p <= 1)
%   Type 2 (Singular Pole)    : \int_0^b (1 / x^p) dx      (Converges if p < 1, Diverges if p >= 1)
%
% Formulas:
%   Type 1 Convergent : I = a^(1 - p) / (p - 1)
%   Type 2 Convergent : I = b^(1 - p) / (1 - p)
%
% Syntax:
%   report = evaluate_p_integral_convergence(p_power, integral_kind, boundary_val)
%
% Inputs:
%   p_power       - Exponent power parameter p (default = 1.5)
%   integral_kind - Integer 1 (for [a, Inf)) or 2 (for [0, b]) (default = 1)
%   boundary_val  - Boundary coordinate a > 0 for Type 1, or b > 0 for Type 2 (default = 1.0)
%
% Outputs:
%   report        - Struct containing convergence classification, closed-form limit, and numerical check

arguments
    p_power       (1,1) double = 1.5
    integral_kind (1,1) double {mustBeMember(integral_kind, [1, 2])} = 1
    boundary_val  (1,1) double {mustBePositive} = 1.0
end

if integral_kind == 1
    % Type 1: \int_a^\infty x^(-p) dx
    model_name = sprintf('Type 1: \\int_{%.2f}^\\infty (1 / x^{%.2f}) dx', boundary_val, p_power);
    is_convergent = (p_power > 1.0);

    if is_convergent
        val_analytic = (boundary_val^(1.0 - p_power)) / (p_power - 1.0);
        f_handle = @(x) x.^(-p_power);
        val_numerical = integral(f_handle, boundary_val, Inf, 'RelTol', 1e-12, 'AbsTol', 1e-14);
    else
        val_analytic = Inf;
        val_numerical = Inf;
    end

else
    % Type 2: \int_0^b x^(-p) dx (Singularity at origin x = 0)
    model_name = sprintf('Type 2: \\int_0^{%.2f} (1 / x^{%.2f}) dx', boundary_val, p_power);
    is_convergent = (p_power < 1.0);

    if is_convergent
        val_analytic = (boundary_val^(1.0 - p_power)) / (1.0 - p_power);
        f_handle = @(x) x.^(-p_power);
        val_numerical = integral(f_handle, 0, boundary_val, 'RelTol', 1e-12, 'AbsTol', 1e-14);
    else
        val_analytic = Inf;
        val_numerical = Inf;
    end
end

if is_convergent
    discrepancy = abs(val_analytic - val_numerical);
else
    discrepancy = NaN;
end

% Assemble structured report
report.integral_type = integral_kind;
report.power_p = p_power;
report.boundary = boundary_val;
report.formulation_name = model_name;
report.is_convergent = is_convergent;
report.closed_form_value = val_analytic;
report.numerical_quadrature_value = val_numerical;
report.discrepancy = discrepancy;

% Console report
fprintf('=== Canonical p-Power Improper Integral Analysis ===\n');
fprintf('Improper Model Formulation : %s\n', model_name);
fprintf('Power Parameter (p)        : %.4f\n', p_power);
fprintf('Convergence Classification : %s\n', mat2str(is_convergent));
if is_convergent
    fprintf('Closed-Form Analytic Value : %.10f\n', val_analytic);
    fprintf('Numerical Quadrature Check : %.10f (Diff: %.2e)\n', val_numerical, discrepancy);
else
    if integral_kind == 1
        fprintf('Divergence Reason          : p = %.2f <= 1.0 (Asymptotic decay too slow)\n', p_power);
    else
        fprintf('Divergence Reason          : p = %.2f >= 1.0 (Pole singularity at origin too steep)\n', p_power);
    end
end
end