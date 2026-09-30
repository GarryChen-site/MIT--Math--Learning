function sol = solve_parabolic_mirror_sensitivity(theta_measured_rad, delta_theta_rad)
% SOLVE_PARABOLIC_MIRROR_SENSITIVITY Evaluates ray-tracing sensitivity for parabolic
% telescope mirror y = x^2 - 1/4 with focal point at origin (PS3 Part II Problem 4).
% Formulas:
%   Ray hits at (a, b) where b = a^2 - 1/4
%   Ray to focus has angle theta: tan(theta) = b / a = a - 1 / (4*a)
%   Derivative: d(tan theta)/da = sec^2(theta) * (dtheta/da) = 1 + 1 / (4*a^2)
%   Sensitivity: dtheta/da = cos^2(theta) * [ 1 + 1 / (4*a^2) ]
%
% Syntax:
%   sol = solve_parabolic_mirror_sensitivity(theta_measured_rad, delta_theta_rad)

arguments
    theta_measured_rad (1,1) double = -pi / 6.0
    delta_theta_rad    (1,1) double {mustBePositive} = 1e-3
end

syms a_sym real
% tan(theta) = a - 1/(4*a) => 4*a^2 - 4*tan(theta)*a - 1 = 0
tan_val = tan(theta_measured_rad);
a_solutions = double(solve(a_sym - 1.0 / (4.0 * a_sym) == tan_val, a_sym));

% Select physical branch (for theta = -pi/6, a > 0 is chosen)
a_operating = a_solutions(a_solutions > 0);
if isempty(a_operating)
    a_operating = a_solutions(1);
else
    a_operating = a_operating(1);
end

b_operating = a_operating^2 - 0.25;

% Sensitivity derivation
sec2_theta = 1.0 + tan_val^2; % 1 / cos^2(theta)
dtan_da = 1.0 + 1.0 / (4.0 * (a_operating^2));
dtheta_da = dtan_da / sec2_theta;
da_dtheta = 1.0 / dtheta_da;

% Position uncertainty
delta_a = abs(da_dtheta) * delta_theta_rad;

% Assemble structured output
sol.measured_angle_rad = theta_measured_rad;
sol.measured_angle_deg = rad2deg(theta_measured_rad);
sol.angle_uncertainty_rad = delta_theta_rad;
sol.mirror_hit_coordinate_a = a_operating;
sol.mirror_hit_coordinate_b = b_operating;
sol.angular_sensitivity_dtheta_da = dtheta_da;
sol.inverse_sensitivity_da_dtheta = da_dtheta;
sol.estimated_position_error_a = delta_a;

% Console report
fprintf('=== PS3 Problem 4: Parabolic Mirror Observation Sensitivity ===\n');
fprintf('Measured Star Angle (theta) : %.4f rad (%.2f deg)\n', ...
    theta_measured_rad, sol.measured_angle_deg);
fprintf('Angular Measurement Noise   : +/-%.2e rad\n', delta_theta_rad);
fprintf('Mirror Impingement Point (a): %.6f (Height b = %.6f)\n', a_operating, b_operating);
fprintf('Sensitivity Factor dtheta/da: %.6f rad/unit\n', dtheta_da);
fprintf('Inverse Sensitivity da/dtheta: %.6f unit/rad\n', da_dtheta);
fprintf('Estimated Ray Position Error: +/-%.6e units\n', delta_a);
end