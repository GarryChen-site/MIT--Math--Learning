function sol = solve_conical_tank_filling(tank_radius_top, tank_height_total, inflow_rate, eval_height)
% SOLVE_CONICAL_TANK_FILLING Calculates instantaneous rate of water level rise dh/dt
% and cross-sectional parameters for an inverted conical container.
% Volume Model : V(h) = (1/3) * pi * ( (R/H)*h )^2 * h = (pi/3) * (R/H)^2 * h^3
% Rate Model   : dV/dt = pi * (R/H)^2 * h^2 * (dh/dt) => dh/dt = (dV/dt) / [ pi * (R/H)^2 * h^2 ]
%
% Syntax:
%   sol = solve_conical_tank_filling(tank_radius_top, tank_height_total, inflow_rate, eval_height)
%
% Inputs:
%   tank_radius_top   - Radius of conical tank top R (ft or m, default = 4.0)
%   tank_height_total - Total height of tank H (ft or m, default = 10.0)
%   inflow_rate       - Volumetric inflow rate dV/dt (ft^3/min or m^3/min, default = 2.0)
%   eval_height       - Instantaneous water depth h (must be <= H, default = 5.0)
%
% Outputs:
%   sol               - Struct containing instantaneous radius, volume, area, and rise rate

arguments
    tank_radius_top   (1,1) double {mustBePositive} = 4.0
    tank_height_total (1,1) double {mustBePositive} = 10.0
    inflow_rate       (1,1) double {mustBePositive} = 2.0
    eval_height       (1,1) double {mustBePositive} = 5.0
end

if eval_height > tank_height_total
    error('Evaluated water height h cannot exceed total tank height H.');
end

% 1. Similar triangles ratio: r(h) = (R / H) * h
shape_ratio = tank_radius_top / tank_height_total;
current_radius = shape_ratio * eval_height;

% 2. Cross-sectional liquid surface area at height h: A(h) = pi * r^2
surface_area = pi * (current_radius^2);

% 3. Current accumulated liquid volume: V(h) = (1/3) * A(h) * h
current_volume = (1.0 / 3.0) * surface_area * eval_height;
total_tank_volume = (1.0 / 3.0) * pi * (tank_radius_top^2) * tank_height_total;

% 4. Rate of water level rise via related rates: dh/dt = (dV/dt) / A(h)
dh_dt = inflow_rate / surface_area;

% Assemble solution structure
sol.tank_top_radius = tank_radius_top;
sol.tank_total_height = tank_height_total;
sol.geometric_ratio_r_over_h = shape_ratio;
sol.inflow_rate_dV_dt = inflow_rate;
sol.evaluated_water_height = eval_height;
sol.current_water_radius = current_radius;
sol.current_surface_area = surface_area;
sol.current_volume = current_volume;
sol.total_capacity = total_tank_volume;
sol.rate_of_level_rise_dh_dt = dh_dt;

% Console report
fprintf('=== Conical Tank Filling Kinematics ===\n');
fprintf('Tank Dimensions             : Top Radius R = %.2f, Height H = %.2f (r/h = %.2f)\n', ...
    tank_radius_top, tank_height_total, shape_ratio);
fprintf('Volumetric Inflow Rate dV/dt: %.4f per unit time\n', inflow_rate);
fprintf('Evaluated Water Height (h)  : %.2f\n', eval_height);
fprintf('Current Water Radius (r)    : %.4f\n', current_radius);
fprintf('Instantaneous Cross Surface : %.4f\n', surface_area);
fprintf('Accumulated Water Volume    : %.4f / %.4f (Capacity Fill: %.1f%%)\n', ...
    current_volume, total_tank_volume, (current_volume / total_tank_volume) * 100);
fprintf('Instantaneous Rise Rate dh/dt: %.6f per unit time (Theoretical: %.6f)\n', ...
    dh_dt, inflow_rate / (pi * (shape_ratio * eval_height)^2));
end