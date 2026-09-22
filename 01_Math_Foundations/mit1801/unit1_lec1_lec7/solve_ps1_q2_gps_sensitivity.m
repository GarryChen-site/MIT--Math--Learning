function report = solve_ps1_q2_gps_sensitivity(s_alt, h0_val, delta_h_list)
% SOLVE_PS1_Q2_GPS_SENSITIVITY Computes finite difference quotients Delta L / Delta h
% and determines error bound constant C such that |Delta L| <= C * |Delta h|.
%
% Syntax:
%   report = solve_ps1_q2_gps_sensitivity(s_alt, h0_val, delta_h_list)
%
% Inputs:
%   s_alt        - Satellite altitude (default = 20000 km)
%   h0_val       - Operating slant distance (e.g., 25000 or 20001)
%   delta_h_list - Array of perturbation steps (default = [1.0, 0.1, 0.01])
%
% Outputs:
%   report       - Struct containing theoretical sensitivity and empirical bounds

arguments
    s_alt        (1,1) double {mustBePositive} = 20000.0
    h0_val       (1,1) double {mustBePositive} = 25000.0
    delta_h_list (1,:) double {mustBePositive} = [1.0, 0.1, 0.01]
end

L0 = sqrt(h0_val^2 - s_alt^2);
analytic_sensitivity = h0_val / L0; % dL/dh

n_cases = length(delta_h_list);
diff_quotients_plus = zeros(n_cases, 1);
diff_quotients_minus = zeros(n_cases, 1);

for k = 1:n_cases
    dh = delta_h_list(k);
    L_plus = sqrt((h0_val + dh)^2 - s_alt^2);
    L_minus = sqrt((h0_val - dh)^2 - s_alt^2);

    diff_quotients_plus(k) = (L_plus - L0) / dh;
    diff_quotients_minus(k) = (L0 - L_minus) / dh;
end

max_empirical_ratio = max([diff_quotients_plus; diff_quotients_minus]);
% Suggested simplest round number bounding constant C
suggested_C = ceil(max_empirical_ratio * 10) / 10;
if suggested_C > 10
    suggested_C = ceil(max_empirical_ratio);
end

report.s_alt = s_alt;
report.h0 = h0_val;
report.L0 = L0;
report.analytic_derivative = analytic_sensitivity;
report.delta_h_tested = delta_h_list;
report.ratios_positive_dh = diff_quotients_plus;
report.ratios_negative_dh = diff_quotients_minus;
report.recommended_bound_C = suggested_C;

% Console report
fprintf('=== Problem 2: GPS Slant Sensitivity Analysis ===\n');
fprintf('Altitude s = %.1f km | Operating Slant Range h0 = %.1f km\n', s_alt, h0_val);
fprintf('Nominal Ground Distance L0 : %.4f km\n', L0);
fprintf('Analytic Derivative dL/dh  : %.4f\n', analytic_sensitivity);
for k = 1:n_cases
    fprintf('  Delta h = %6.2f -> (+dh): %.5f | (-dh): %.5f\n', ...
        delta_h_list(k), diff_quotients_plus(k), diff_quotients_minus(k));
end
fprintf('Simplest Round Constant C  : |Delta L| <= %.1f * |Delta h|\n', suggested_C);
end