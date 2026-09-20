function threshold_info = solve_trig_linearization_threshold(max_tol_pct)
% SOLVE_TRIG_LINEARIZATION_THRESHOLD Solves for maximum angle bounds
% where sin(theta) approx theta holds within a relative tolerance percentage.
%
% Syntax:
%   threshold_info = solve_trig_linearization_threshold(max_tol_pct)
%
% Inputs:
%   max_tol_pct - Maximum acceptable relative error percentage (e.g., 1.0 for 1%)
%
% Outputs:
%   threshold_info - Struct containing valid angle limits in radians and degrees

arguments
    max_tol_pct (1,1) double {mustBePositive} = 1.0
end

% Relative error equation: (theta - sin(theta)) / theta = tol
tol_ratio = max_tol_pct / 100.0;
residual_eqn = @(th) (th - sin(th)) ./ th - tol_ratio;

% Solve numerically using fzero with an initial guess
theta_limit_rad = fzero(residual_eqn, [0.01, pi/2]);
theta_limit_deg = rad2deg(theta_limit_rad);

threshold_info.tolerance_pct = max_tol_pct;
threshold_info.max_angle_rad = theta_limit_rad;
threshold_info.max_angle_deg = theta_limit_deg;

% Console report
fprintf('=== Small-Angle Linearization Validity Threshold ===\n');
fprintf('Design Tolerance Limit    : %.2f%%\n', max_tol_pct);
fprintf('Maximum Angle (Radians)   : %.4f rad\n', theta_limit_rad);
fprintf('Maximum Angle (Degrees)   : %.2f deg\n', theta_limit_deg);
fprintf('Rule-of-Thumb Applicability: Linearization holds within [-%0.1f deg, +%0.1f deg]\n', ...
    theta_limit_deg, theta_limit_deg);
end