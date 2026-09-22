function result = solve_ps1_q5_tank_drainage(drain_time_limit, denom_coeff, eval_t)
% SOLVE_PS1_Q5_TANK_DRAINAGE Computes average discharge rate over initial period
% and instantaneous discharge rate at time t for volume model: V(t) = (T - t)^2 / denom_coeff.
%
% Syntax:
%   result = solve_ps1_q5_tank_drainage(drain_time_limit, denom_coeff, eval_t)
%
% Inputs:
%   drain_time_limit - Complete dry-out time parameter T (minutes, default = 10)
%   denom_coeff      - Denominator scalar coefficient (default = 5)
%   eval_t           - Specific time evaluation point (minutes, default = 5)
%
% Outputs:
%   result           - Struct containing initial volume, average rate, and instantaneous flow rate

arguments
    drain_time_limit (1,1) double {mustBePositive} = 10.0
    denom_coeff      (1,1) double {mustBePositive} = 5.0
    eval_t           (1,1) double {mustBeNonnegative} = 5.0
end

V = @(t) ((drain_time_limit - t).^2) ./ denom_coeff;
dV_dt = @(t) -2 * (drain_time_limit - t) ./ denom_coeff;

V0 = V(0);
V_at_eval = V(eval_t);

% Average drainage rate: [V(eval_t) - V(0)] / eval_t
avg_rate_of_change = (V_at_eval - V0) / eval_t;
avg_drain_speed = abs(avg_rate_of_change);

% Instantaneous rate at eval_t
inst_rate_of_change = dV_dt(eval_t);
inst_drain_speed = abs(inst_rate_of_change);

result.initial_volume = V0;
result.volume_at_t = V_at_eval;
result.average_rate_of_change = avg_rate_of_change;
result.average_drain_rate = avg_drain_speed;
result.instantaneous_rate_of_change = inst_rate_of_change;
result.instantaneous_drain_rate = inst_drain_speed;

% Console report
fprintf('=== Problem 5: Water Cooler Drainage Kinematics ===\n');
fprintf('Initial Volume at t=0      : %.2f L\n', V0);
fprintf('Volume at t = %.1f min       : %.2f L\n', eval_t, V_at_eval);
fprintf('Average Drainage Rate (0-5m): %.2f L/min (Rate of Change: %.2f L/min)\n', ...
    avg_drain_speed, avg_rate_of_change);
fprintf('Instantaneous Rate at t=%.1fm : %.2f L/min (Rate of Change: %.2f L/min)\n', ...
    eval_t, inst_drain_speed, inst_rate_of_change);
end