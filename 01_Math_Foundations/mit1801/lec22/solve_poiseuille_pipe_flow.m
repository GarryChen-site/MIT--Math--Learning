function report = solve_poiseuille_pipe_flow(pipe_radius_R, visc_coeff_c, delta_R_percent)
% SOLVE_POISEUILLE_PIPE_FLOW Analyzes Poiseuille fluid flow through a cylindrical pipe:
% Velocity Field : v(r) = c * (R^2 - r^2)
% Total Flow     : Q = \int_0^R v(r) * (2*pi*r dr) = (pi / 2) * c * R^4
% Mean Velocity  : v_avg = Q / (pi * R^2) = 0.5 * v_max
% Sensitivity    : dQ / dR = 2 * pi * c * R^3 => (dQ / Q) = 4 * (dR / R)
%
% Syntax:
%   report = solve_poiseuille_pipe_flow(pipe_radius_R, visc_coeff_c, delta_R_percent)
%
% Inputs:
%   pipe_radius_R   - Pipe inner radius R in meters (default = 0.05 m)
%   visc_coeff_c    - Flow velocity factor c in 1/(m*s) (default = 100.0)
%   delta_R_percent - Hypothetical percentage expansion in radius (default = 5.0 for +5%)
%
% Outputs:
%   report          - Struct containing maximum speed, flow rate Q, and sensitivity scaling

arguments
    pipe_radius_R   (1,1) double {mustBePositive} = 0.05
    visc_coeff_c    (1,1) double {mustBePositive} = 100.0
    delta_R_percent (1,1) double = 5.0
end

% 1. Kinematic Velocity Metrics
v_max = visc_coeff_c * (pipe_radius_R^2); % At r = 0
cross_section_area = pi * (pipe_radius_R^2);

% 2. Volumetric Flow Rate Q via Annular Integration
% Q = 2*pi*c * \int_0^R (R^2*r - r^3) dr = (pi / 2) * c * R^4
flow_rate_Q = (pi / 2.0) * visc_coeff_c * (pipe_radius_R^4);
v_mean = flow_rate_Q / cross_section_area; % Exactly 0.5 * v_max

% 3. Sensitivity Analysis (Poiseuille R^4 Law):
% Relative flow amplification: Q_new / Q_old = (1 + delta_R/100)^4
fractional_r_increase = delta_R_percent / 100.0;
flow_amplification_exact = (1.0 + fractional_r_increase)^4 - 1.0;
% Linear differential approximation: dQ/Q approx 4 * (dR/R)
flow_amplification_linear = 4.0 * fractional_r_increase;

% Assemble structured report
report.pipe_radius_m = pipe_radius_R;
report.cross_section_area_m2 = cross_section_area;
report.flow_constant_c = visc_coeff_c;
report.max_center_velocity_mps = v_max;
report.mean_velocity_mps = v_mean;
report.volumetric_flow_rate_m3_s = flow_rate_Q;
report.flow_liters_per_second = flow_rate_Q * 1000.0;
report.radius_perturbation_pct = delta_R_percent;
report.exact_flow_increase_pct = flow_amplification_exact * 100.0;
report.linear_approx_increase_pct = flow_amplification_linear * 100.0;

% Console report
fprintf('=== Poiseuille Laminar Pipe Flow Kinematics ===\n');
fprintf('Pipe Geometry Radius (R)    : %.4f m (Area = %.4e m^2)\n', ...
    pipe_radius_R, cross_section_area);
fprintf('Centerline Peak Speed (v_max): %.4f m/s\n', v_max);
fprintf('Average Flow Speed (v_mean) : %.4f m/s (Theory: 0.5 * v_max)\n', v_mean);
fprintf('Volumetric Flow Rate (Q)    : %.6e m^3/s (%.3f L/s)\n', ...
    flow_rate_Q, report.flow_liters_per_second);
fprintf('--- Poiseuille R^4 Sensitivity Law ---\n');
fprintf('Radius Expansion Input      : %+.2f%%\n', delta_R_percent);
fprintf('Exact Resulting Flow Gain   : %+.2f%%\n', report.exact_flow_increase_pct);
fprintf('1st-Order Sensitivity Bound : %+.2f%% (Rule: 4 * dR/R)\n', ...
    report.linear_approx_increase_pct);
end