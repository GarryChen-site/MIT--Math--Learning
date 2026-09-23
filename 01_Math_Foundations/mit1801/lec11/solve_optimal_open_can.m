function sol = solve_optimal_open_can(target_volume)
% SOLVE_OPTIMAL_OPEN_CAN Computes optimal geometric proportions and minimum surface area
% for an open-topped cylindrical container of fixed volume V using constraint substitution.
% Constraint: V = pi * r^2 * h = const
% Objective : S(r) = pi * r^2 + 2*V/r
%
% Syntax:
%   sol = solve_optimal_open_can(target_volume)
%
% Inputs:
%   target_volume - Target enclosed volume V (e.g., in cm^3 or m^3, default = 1000.0)
%
% Outputs:
%   sol           - Struct containing optimal radius, height, minimum area, and aspect ratio

arguments
    target_volume (1,1) double {mustBePositive} = 1000.0
end

syms r V positive
% 1. Eliminate height h via volume constraint: h = V / (pi * r^2)
h_constraint = V / (pi * r^2);
S_sym = pi * r^2 + 2 * pi * r * h_constraint; % S(r) = pi*r^2 + 2*V/r

% 2. First and second derivatives
dS_dr = diff(S_sym, r);
d2S_dr2 = diff(dS_dr, r);

% 3. Critical point: dS/dr = 0 => 2*pi*r - 2*V/r^2 = 0 => r* = (V/pi)^(1/3)
r_crit_sym = solve(dS_dr == 0, r);
r_opt_val = double(subs(r_crit_sym, V, target_volume));

% 4. Optimal height and minimum area
h_opt_val = double(subs(h_constraint, [r, V], [r_opt_val, target_volume]));
S_min_val = double(subs(S_sym, [r, V], [r_opt_val, target_volume]));
d2S_val   = double(subs(d2S_dr2, [r, V], [r_opt_val, target_volume]));

% Populate solution structure
sol.volume = target_volume;
sol.optimal_radius = r_opt_val;
sol.optimal_height = h_opt_val;
sol.aspect_ratio_h_over_r = h_opt_val / r_opt_val;
sol.minimum_surface_area = S_min_val;
sol.second_derivative_value = d2S_val;
sol.is_convex_minimum = (d2S_val > 0);

% Console report
fprintf('=== Open Cylindrical Can Surface Area Optimization ===\n');
fprintf('Fixed Volume Capacity (V) : %.2f\n', target_volume);
fprintf('Optimal Base Radius (r*)  : %.4f\n', r_opt_val);
fprintf('Optimal Can Height (h*)   : %.4f\n', h_opt_val);
fprintf('Optimal Aspect Ratio (h/r): %.4f (Theoretical: 1.0000)\n', sol.aspect_ratio_h_over_r);
fprintf('Minimum Surface Area (S*) : %.4f\n', S_min_val);
fprintf('Second Derivative d2S/dr2 : %.4f (Strictly Convex Min: %s)\n', ...
    d2S_val, mat2str(sol.is_convex_minimum));
end