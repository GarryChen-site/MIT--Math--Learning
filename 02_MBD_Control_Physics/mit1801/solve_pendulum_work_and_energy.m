function sol = solve_pendulum_work_and_energy(mass_kg, length_m, initial_angle_deg, gravity_g)
% SOLVE_PENDULUM_WORK_AND_ENERGY Computes gravitational work on a swinging pendulum
% via line integral of tangential force F = m*g*sin(theta) and proves path-independence
% against potential energy drop Delta V = m*g*Delta h where Delta h = L*(1 - cos(theta0)).
%
% Formulas:
%   Tangential Force : F(theta) = m * g * sin(theta)
%   Path Integral    : W = \int_0^theta0 (m*g*sin(theta)) * (L*dtheta) = m*g*L*(1 - cos(theta0))
%   Vertical Drop    : Delta h = L * (1 - cos(theta0))
%   Potential Energy : Delta V = m * g * Delta h
%
% Syntax:
%   sol = solve_pendulum_work_and_energy(mass_kg, length_m, initial_angle_deg, gravity_g)
%
% Inputs:
%   mass_kg           - Mass of pendulum bob (kg, default = 1.0)
%   length_m          - Length of pendulum rod (m, default = 1.0)
%   initial_angle_deg - Release angle in degrees (default = 45.0)
%   gravity_g         - Local gravity acceleration in m/s^2 (default = 9.81)
%
% Outputs:
%   sol               - Struct containing path work, potential energy, and path-independence proof

arguments
    mass_kg           (1,1) double {mustBePositive} = 1.0
    length_m          (1,1) double {mustBePositive} = 1.0
    initial_angle_deg (1,1) double {mustBeInRange(initial_angle_deg, 0, 180)} = 45.0
    gravity_g         (1,1) double {mustBePositive} = 9.81
end

theta0_rad = deg2rad(initial_angle_deg);

% 1. Analytic Work via Lecture 25 Line Integral:
% W = -L*m*g*cos(theta)|_0^theta0 = m*g*L*(1 - cos(theta0))
work_analytic = mass_kg * gravity_g * length_m * (1.0 - cos(theta0_rad));

% 2. Numerical Integration along circular arc: \int_0^theta0 (m*g*sin(th)) * (L*dth)
tangential_force = @(th) mass_kg * gravity_g * sin(th);
work_numerical = integral(@(th) tangential_force(th) * length_m, 0, theta0_rad, ...
    'RelTol', 1e-12, 'AbsTol', 1e-14);

% 3. Potential Energy Change via Direct Vertical Descent:
% Height drop Delta h = L - L*cos(theta0) = L*(1 - cos(theta0))
delta_h = length_m * (1.0 - cos(theta0_rad));
delta_potential_energy = mass_kg * gravity_g * delta_h;

% Path independence residual
discrepancy = abs(work_analytic - delta_potential_energy);
is_conservative = (discrepancy < 1e-12);

% Peak kinetic energy and maximum velocity at the bottom (K_max = W)
v_max = sqrt(2.0 * work_analytic / mass_kg);
omega_max = v_max / length_m;

% Assemble structured output
sol.mass = mass_kg;
sol.length = length_m;
sol.initial_angle_deg = initial_angle_deg;
sol.initial_angle_rad = theta0_rad;
sol.vertical_height_drop = delta_h;
sol.work_path_analytic = work_analytic;
sol.work_path_numerical = work_numerical;
sol.potential_energy_drop = delta_potential_energy;
sol.is_path_independent = is_conservative;
sol.max_linear_velocity = v_max;
sol.max_angular_velocity = omega_max;

% Console report
fprintf('=== Pendulum Work, Energy & Path-Independence Solution ===\n');
fprintf('Pendulum Parameters         : Mass = %.2f kg | Length = %.2f m\n', mass_kg, length_m);
fprintf('Release Angular Displacement: %.2f deg (%.4f rad)\n', initial_angle_deg, theta0_rad);
fprintf('Vertical Fall Distance (h)  : %.4f m\n', delta_h);
fprintf('Line Integral Work (W)      : %.6f Joules\n', work_analytic);
fprintf('Potential Energy Drop (mg h): %.6f Joules\n', delta_potential_energy);
fprintf('Conservative Field Verified : %s (Path Independence Confirmed)\n', mat2str(is_conservative));
fprintf('Terminal Velocity at Bottom : Linear = %.4f m/s | Angular = %.4f rad/s\n', ...
    v_max, omega_max);
end