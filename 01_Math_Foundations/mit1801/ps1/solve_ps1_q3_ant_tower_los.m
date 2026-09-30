function sol = solve_ps1_q3_ant_tower_los(H_mound, H_tower, k_shape)
% SOLVE_PS1_Q3_ANT_TOWER_LOS Calculates critical tangency height where an observer
% on a parabolic mound outline y = H_mound - k*x^2 establishes direct line-of-sight
% with the tip of a tower of height H_tower positioned at x = 0.
%
% Syntax:
%   sol = solve_ps1_q3_ant_tower_los(H_mound, H_tower, k_shape)
%
% Inputs:
%   H_mound - Height of the mound at peak (ft, default = 1000)
%   H_tower - Height of the tower on top (ft, default = 100)
%   k_shape - Curvature coefficient of parabola (default = 1.0)
%
% Outputs:
%   sol     - Struct containing tangency coordinates (x0, y0) and LOS slope

arguments
    H_mound (1,1) double {mustBePositive} = 1000.0
    H_tower (1,1) double {mustBePositive} = 100.0
    k_shape (1,1) double {mustBePositive} = 1.0
end

% Tip of the tower coordinate: (0, H_mound + H_tower)
y_tip = H_mound + H_tower;

% Parabola: y(x) = H_mound - k*x^2 => y'(x) = -2*k*x
% Tangent line through (0, y_tip) and (x0, y(x0)):
% Slope equality: (y_tip - y(x0)) / (0 - x0) = -2*k*x0
% => [H_tower + k*x0^2] / (-x0) = -2*k*x0
% => H_tower + k*x0^2 = 2*k*x0^2 => k*x0^2 = H_tower => x0 = sqrt(H_tower / k)
x0_crit = sqrt(H_tower / k_shape);
y0_crit = H_mound - k_shape * x0_crit^2; % Exactly: H_mound - H_tower
slope_crit = -2 * k_shape * x0_crit;

sol.tower_tip_height = y_tip;
sol.critical_x = x0_crit;
sol.critical_height_y = y0_crit;
sol.los_tangent_slope = slope_crit;

% Console report
fprintf('=== Problem 3: Line-of-Sight Visibility Threshold ===\n');
fprintf('Mound Peak Height   : %.1f ft\n', H_mound);
fprintf('Tower Height        : %.1f ft (Tower Tip at y = %.1f ft)\n', H_tower, y_tip);
fprintf('Critical Position   : x0 = %.2f ft\n', x0_crit);
fprintf('Visibility Height y : %.2f ft (Ant starts seeing tower here)\n', y0_crit);
fprintf('LOS Tangent Slope   : %.4f\n', slope_crit);
end