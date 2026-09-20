function sol = solve_gps_sensitivity(s_alt, h_measured, delta_h)
% SOLVE_GPS_SENSITIVITY Analyzes GPS flat-earth model sensitivity:
% Relationship: h^2 = s^2 + L^2
% Sensitivity derivative: dL/dh = h / L
%
% Syntax:
%   sol = solve_gps_sensitivity(s_alt, h_measured, delta_h)
%
% Inputs:
%   s_alt      - Satellite altitude above ground (e.g., in km)
%   h_measured - Measured slant range (must be > s_alt)
%   delta_h    - Slant range measurement uncertainty (same units as h)
%
% Outputs:
%   sol        - Struct containing ground distance, sensitivity, and position uncertainty

arguments
    s_alt      (1,1) double {mustBePositive} = 20000.0
    h_measured (1,1) double {mustBePositive} = 25000.0
    delta_h    (1,1) double {mustBeNonnegative} = 0.010 % 10 meters uncertainty
end

if h_measured <= s_alt
    error('Measured slant range h must be strictly greater than satellite altitude s.');
end

% 1. Compute ground distance L = sqrt(h^2 - s^2)
L_val = sqrt(h_measured^2 - s_alt^2);

% 2. Derivative (Sensitivity factor): dL/dh = h / L
sensitivity_factor = h_measured / L_val;

% 3. Error propagation: Delta L approx |dL/dh| * Delta h
delta_L = sensitivity_factor * delta_h;

% Assemble solution struct
sol.satellite_altitude = s_alt;
sol.slant_range = h_measured;
sol.ground_distance = L_val;
sol.sensitivity_dL_dh = sensitivity_factor;
sol.measurement_error_h = delta_h;
sol.position_uncertainty_L = delta_L;

% Console report
fprintf('=== GPS Sensitivity and Error Propagation Analysis ===\n');
fprintf('Satellite Altitude (s)       : %.3f km\n', s_alt);
fprintf('Measured Slant Range (h)     : %.3f km\n', h_measured);
fprintf('Computed Ground Distance (L) : %.3f km\n', L_val);
fprintf('Sensitivity Factor (dL/dh)   : %.4f\n', sensitivity_factor);
fprintf('Slant Range Error (Delta h)  : %.6f km (%.2f m)\n', delta_h, delta_h * 1000);
fprintf('Ground Uncertainty (Delta L) : %.6f km (%.2f m)\n', delta_L, delta_L * 1000);
end