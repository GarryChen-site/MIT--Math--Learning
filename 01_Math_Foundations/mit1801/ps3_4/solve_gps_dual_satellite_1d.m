function sol = solve_gps_dual_satellite_1d(s_alt, x1, delta_x1, h1, delta_h1, x2, delta_x2, h2, delta_h2)
% SOLVE_GPS_DUAL_SATELLITE_1D Solves user position and propagated uncertainty
% on a 1D flat planet using two hovering satellites at altitude s (PS3 Part II Problem 3).
% Models:
%   Ground Distances  : L1 = sqrt(h1^2 - s^2),  L2 = sqrt(h2^2 - s^2)
%   Sensitivity Gains : dL1/dh1 = h1 / L1,      dL2/dh2 = h2 / L2
%   Error Propagation : Delta L_i = (h_i / L_i) * Delta h_i + Delta x_i
%
% Syntax:
%   sol = solve_gps_dual_satellite_1d(s_alt, x1, delta_x1, h1, delta_h1, x2, delta_x2, h2, delta_h2)

arguments
    s_alt    (1,1) double {mustBePositive} = 20000.0
    x1       (1,1) double = 0.0
    delta_x1 (1,1) double {mustBeNonnegative} = 1e-10
    h1       (1,1) double {mustBePositive} = 21000.0
    delta_h1 (1,1) double {mustBeNonnegative} = 1e-2
    x2       (1,1) double = 41596.88
    delta_x2 (1,1) double {mustBeNonnegative} = 1e-10
    h2       (1,1) double {mustBePositive} = 52000.0
    delta_h2 (1,1) double {mustBeNonnegative} = 1e-2
end

if h1 <= s_alt || h2 <= s_alt
    error('Measured slant ranges h1 and h2 must be strictly greater than altitude s.');
end

% 1. Radial Ground Radii from Sub-Satellite Points
L1 = sqrt(h1^2 - s_alt^2); % approx 6403.12 km
L2 = sqrt(h2^2 - s_alt^2); % approx 48000.00 km

% 2. Sensitivities (Jacobian elements) and Elevation Angles
dL1_dh1 = h1 / L1;
dL2_dh2 = h2 / L2;

theta1_deg = rad2deg(asin(s_alt / h1));
theta2_deg = rad2deg(asin(s_alt / h2));

% 3. Individual Uncertainties
delta_L1_total = dL1_dh1 * delta_h1 + delta_x1;
delta_L2_total = dL2_dh2 * delta_h2 + delta_x2;

% 4. Candidate Locations in 1D Space: x_user = x_i +/- L_i
cands_sat1 = [x1 - L1, x1 + L1];
cands_sat2 = [x2 - L2, x2 + L2];

% Intersect candidate sets to find matching location
tol_match = 1.0; % 1 km resolution
matched_pos = NaN;
for p1 = cands_sat1
    for p2 = cands_sat2
        if abs(p1 - p2) < tol_match
            matched_pos = (p1 + p2) / 2.0;
        end
    end
end

% Optimal Variance Weighted Location Uncertainty (Sensor Fusion)
sigma_fused = sqrt(1.0 / ( (1.0 / delta_L1_total^2) + (1.0 / delta_L2_total^2) ));

% Assemble solution structure
sol.ground_radius_L1 = L1;
sol.ground_radius_L2 = L2;
sol.sensitivity_dL1_dh1 = dL1_dh1;
sol.sensitivity_dL2_dh2 = dL2_dh2;
sol.elevation_angles_deg = [theta1_deg, theta2_deg];
sol.position_uncertainty_L1 = delta_L1_total;
sol.position_uncertainty_L2 = delta_L2_total;
sol.resolved_user_position = matched_pos;
sol.fused_position_accuracy = sigma_fused;

% Console report
fprintf('=== PS3 Problem 3: 1D Dual-Satellite GDOP Positioning ===\n');
fprintf('Sat 1: Sub-point = %.2f km, Slant h1 = %.2f km -> L1 = %.2f km\n', x1, h1, L1);
fprintf('  Sensitivity dL1/dh1 : %.4f (Elev Angle: %.2f deg) | Uncertainty: +/-%.4e km\n', ...
    dL1_dh1, theta1_deg, delta_L1_total);
fprintf('Sat 2: Sub-point = %.2f km, Slant h2 = %.2f km -> L2 = %.2f km\n', x2, h2, L2);
fprintf('  Sensitivity dL2/dh2 : %.4f (Elev Angle: %.2f deg) | Uncertainty: +/-%.4e km\n', ...
    dL2_dh2, theta2_deg, delta_L2_total);
fprintf('Resolved User Position : x_user = %.4f km\n', matched_pos);
fprintf('Fused Measurement Acc  : +/-%.4e km (+/-%.2f meters)\n', ...
    sigma_fused, sigma_fused * 1000.0);
end