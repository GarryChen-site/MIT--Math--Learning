function sol = solve_freefall_kinematics(h0, c)
% SOLVE_FREEFALL_KINEMATICS Computes impact timing, average velocity, and
% instantaneous impact velocity for free-fall model: y(t) = h0 - c*t^2.
%
% Syntax:
%   sol = solve_freefall_kinematics(h0, c)
%
% Inputs:
%   h0 - Initial height (ft or m, default = 400.0)
%   c  - Acceleration coefficient (1/2 * g, default = 16.0 for imperial)
%
% Outputs:
%   sol - Struct containing impact time, average velocity, and impact velocity

arguments
    h0 (1,1) double {mustBePositive} = 400.0
    c  (1,1) double {mustBePositive} = 16.0
end

% 1. Ground impact condition: y(t) = 0 => h0 - c * t^2 = 0
t_impact = sqrt(h0 / c);

% 2. Average velocity across the entire motion: Delta y / Delta t
v_average = (0 - h0) / t_impact;

% 3. Instantaneous velocity: v(t) = dy/dt = -2 * c * t
v_impact = -2 * c * t_impact;

% Structure assembly
sol.initial_height = h0;
sol.acceleration_coeff = c;
sol.impact_time = t_impact;
sol.average_velocity = v_average;
sol.instantaneous_impact_velocity = v_impact;
sol.impact_speed_mph = abs(v_impact) * 3600 / 5280;

% Console report
fprintf('=== Free-Fall Kinematics Solution ===\n');
fprintf('Initial Height    : %.2f ft\n', h0);
fprintf('Impact Time       : %.4f s\n', t_impact);
fprintf('Average Velocity  : %.2f ft/s\n', v_average);
fprintf('Instantaneous Vel : %.2f ft/s (%.2f mph)\n', v_impact, sol.impact_speed_mph);
end