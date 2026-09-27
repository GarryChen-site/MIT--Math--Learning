function sol = solve_ftc_kinematics(v_handle, t_interval, initial_pos_x0)
% SOLVE_FTC_KINEMATICS Solves the First Fundamental Theorem of Calculus for kinematics:
% Computes net displacement, odometer total distance, average velocity, and average speed.
%
% Formulas:
%   Net Displacement : Delta x = \int_t0^tf v(t) dt = x(tf) - x(t0)
%   Total Distance   : S = \int_t0^tf |v(t)| dt
%   Final Position   : x(tf) = x0 + Delta x
%
% Syntax:
%   sol = solve_ftc_kinematics(v_handle, t_interval, initial_pos_x0)
%
% Inputs:
%   v_handle        - Function handle for velocity v(t) (default: @(t) sin(t))
%   t_interval      - Two-element vector [t0, tf] with t0 < tf (default: [0, 2*pi])
%   initial_pos_x0  - Initial position x(t0) (default = 0.0)
%
% Outputs:
%   sol             - Struct containing displacement, distance, position, and averages

arguments
    v_handle       function_handle = @(t) sin(t)
    t_interval     (1,2) double = [0.0, 2*pi]
    initial_pos_x0 (1,1) double = 0.0
end

t0 = t_interval(1);
tf = t_interval(2);

if t0 >= tf
    error('Initial time t0 must be strictly less than final time tf.');
end

duration = tf - t0;

% 1. Net displacement via FTC 1: \int_t0^tf v(t) dt
net_disp = integral(v_handle, t0, tf, 'RelTol', 1e-12, 'AbsTol', 1e-14);

% 2. Total distance traveled: \int_t0^tf |v(t)| dt
v_abs = @(t) abs(v_handle(t));
total_dist = integral(v_abs, t0, tf, 'RelTol', 1e-12, 'AbsTol', 1e-14);

% 3. Final state position
final_pos = initial_pos_x0 + net_disp;

% 4. Kinematic averages
avg_velocity = net_disp / duration;
avg_speed = total_dist / duration;

% Assemble solution structure
sol.time_interval = [t0, tf];
sol.duration_seconds = duration;
sol.initial_position = initial_pos_x0;
sol.net_displacement = net_disp;
sol.final_position = final_pos;
sol.total_distance_traveled = total_dist;
sol.average_velocity = avg_velocity;
sol.average_speed = avg_speed;
sol.cancellation_ratio = (total_dist - abs(net_disp)) / total_dist;

% Console report
fprintf('=== FTC 1 Kinematics & Displacement Solver ===\n');
fprintf('Time Window [t0, tf]         : [%.4f, %.4f] s (Duration: %.4f s)\n', t0, tf, duration);
fprintf('Initial State Position x(t0) : %.4f\n', initial_pos_x0);
fprintf('Net Displacement (FTC 1)     : %.8f\n', net_disp);
fprintf('Final State Position x(tf)   : %.8f\n', final_pos);
fprintf('Odometer Total Distance (S)  : %.8f\n', total_dist);
fprintf('Average Velocity (dx/dt_avg) : %.6f per second\n', avg_velocity);
fprintf('Average Speed (|v|_avg)      : %.6f per second\n', avg_speed);
fprintf('Direction Cancellation Ratio : %.2f%%\n', sol.cancellation_ratio * 100);
end