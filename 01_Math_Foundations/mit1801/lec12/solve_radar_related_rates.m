function sol = solve_radar_related_rates(y_offset, D_measured, dD_dt_measured, speed_limit_mph)
% SOLVE_RADAR_RELATED_RATES Computes actual vehicle speed along a straight road
% from line-of-sight radar range and range rate using Pythagorean related rates.
% Geometry : y_offset^2 + x^2 = D^2
% Rates    : 2*x*(dx/dt) = 2*D*(dD/dt) => dx/dt = (D / x) * (dD/dt)
%
% Syntax:
%   sol = solve_radar_related_rates(y_offset, D_measured, dD_dt_measured, speed_limit_mph)
%
% Inputs:
%   y_offset         - Perpendicular distance from radar to road (ft or m, default = 30.0)
%   D_measured       - Measured slant range to target vehicle (must be > y_offset, default = 50.0)
%   dD_dt_measured   - Measured radial range rate (negative if approaching, default = -80.0)
%   speed_limit_mph  - Regulatory speed limit in miles per hour (default = 65.0)
%
% Outputs:
%   sol              - Struct containing horizontal coordinate, true velocity, and speeding verdict

arguments
    y_offset        (1,1) double {mustBePositive} = 30.0
    D_measured      (1,1) double {mustBePositive} = 50.0
    dD_dt_measured  (1,1) double = -80.0
    speed_limit_mph (1,1) double {mustBePositive} = 65.0
end

if D_measured <= y_offset
    error('Measured slant range D must be strictly greater than perpendicular road offset.');
end

% 1. Horizontal position along road: x = sqrt(D^2 - y^2)
x_pos = sqrt(D_measured^2 - y_offset^2);

% 2. True car velocity along the road via related rates
% Derivative of x^2 + y0^2 = D^2 with y0 constant: 2*x*(dx/dt) = 2*D*(dD/dt)
dx_dt = (D_measured / x_pos) * dD_dt_measured;

% 3. Unit conversion and speed limit compliance
speed_true_fps = abs(dx_dt);
speed_true_mph = speed_true_fps * 3600.0 / 5280.0;
speed_limit_fps = speed_limit_mph * 5280.0 / 3600.0;
is_speeding = (speed_true_mph > speed_limit_mph);

% Assemble solution structure
sol.perpendicular_offset_ft = y_offset;
sol.slant_range_D_ft = D_measured;
sol.range_rate_dD_dt_fps = dD_dt_measured;
sol.horizontal_distance_x_ft = x_pos;
sol.projection_factor = D_measured / x_pos; % 1 / cos(theta)
sol.actual_velocity_dx_dt_fps = dx_dt;
sol.actual_speed_mph = speed_true_mph;
sol.speed_limit_mph = speed_limit_mph;
sol.speed_limit_fps = speed_limit_fps;
sol.is_speeding = is_speeding;

% Console report
fprintf('=== Police Radar Tracking & Speeding Detection ===\n');
fprintf('Perpendicular Road Offset  : %.2f ft\n', y_offset);
fprintf('Measured Slant Distance (D): %.2f ft\n', D_measured);
fprintf('Measured Range Rate (dD/dt): %.2f ft/s\n', dD_dt_measured);
fprintf('Vehicle Road Position (x)  : %.2f ft\n', x_pos);
fprintf('Geometric Scaling (D/x)    : %.4f (1 / cos(\\theta))\n', sol.projection_factor);
fprintf('Calculated Actual Velocity : %.2f ft/s (%.2f mph)\n', dx_dt, speed_true_mph);
fprintf('Legal Speed Limit          : %.2f mph (%.2f ft/s)\n', speed_limit_mph, speed_limit_fps);
fprintf('Speeding Infraction Status : %s (Excess: %+.2f mph)\n', ...
    mat2str(is_speeding), speed_true_mph - speed_limit_mph);
end