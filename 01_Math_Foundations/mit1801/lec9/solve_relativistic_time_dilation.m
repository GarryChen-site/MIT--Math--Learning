function report = solve_relativistic_time_dilation(v_orbital, c_light, nominal_time_T)
% SOLVE_RELATIVISTIC_TIME_DILATION Computes exact and approximated relativistic
% time dilation for orbiting satellites (GPS calibration model).
% Formula: T' = T / sqrt(1 - v^2/c^2) approx T * [ 1 + 0.5*(v/c)^2 + (3/8)*(v/c)^4 ]
%
% Syntax:
%   report = solve_relativistic_time_dilation(v_orbital, c_light, nominal_time_T)
%
% Inputs:
%   v_orbital      - Satellite velocity (km/s, default = 4.0 km/s from lecture)
%   c_light        - Speed of light (km/s, default = 3.0e5 km/s)
%   nominal_time_T - Onboard clock baseline duration (seconds, default = 86400 s, 1 day)
%
% Outputs:
%   report         - Struct containing exact dilation, linear offset, and higher-order magnitude

arguments
    v_orbital      (1,1) double {mustBePositive} = 4.0
    c_light        (1,1) double {mustBePositive} = 3.0e5
    nominal_time_T (1,1) double {mustBePositive} = 86400.0 % 24 hours in seconds
end

beta = v_orbital / c_light;
beta_sq = beta^2;

% 1. Exact Special Relativity Dilation
gamma = 1.0 / sqrt(1.0 - beta_sq);
T_prime_exact = nominal_time_T * gamma;
exact_time_shift = T_prime_exact - nominal_time_T;

% 2. 1st-Order Linear Approximation: gamma approx 1 + 0.5 * beta^2
linear_shift = nominal_time_T * (0.5 * beta_sq);

% 3. 2nd-Order Quadratic Term: (3/8) * beta^4
quadratic_term_fraction = (3.0 / 8.0) * (beta_sq^2);
quadratic_shift = nominal_time_T * quadratic_term_fraction;

% Assemble report structure
report.velocity_km_s = v_orbital;
report.speed_of_light = c_light;
report.beta_ratio = beta;
report.beta_squared = beta_sq;
report.exact_clock_shift_sec = exact_time_shift;
report.linear_approx_shift_sec = linear_shift;
report.quadratic_offset_sec = quadratic_shift;
report.is_quadratic_negligible = (quadratic_shift < 1e-15); % Below typical atomic clock threshold

% Console report
fprintf('=== GPS Satellite Relativistic Time Dilation Analysis ===\n');
fprintf('Satellite Velocity (v)     : %.2f km/s\n', v_orbital);
fprintf('Speed of Light (c)         : %.2e km/s\n', c_light);
fprintf('Dimensionless (v/c)^2      : %.4e\n', beta_sq);
fprintf('Nominal Duration (T)       : %.1f s (24 hours)\n', nominal_time_T);
fprintf('Exact Clock Dilation Offset: %.8e s (%.3f microseconds/day)\n', ...
    exact_time_shift, exact_time_shift * 1e6);
fprintf('1st-Order Linear Correction: %.8e s (%.3f microseconds/day)\n', ...
    linear_shift, linear_shift * 1e6);
fprintf('2nd-Order Quadratic Term   : %.8e s\n', quadratic_shift);
fprintf('Quadratic Term Negligible? : %s (Requires atomic precision < 1e-15 s)\n', ...
    mat2str(report.is_quadratic_negligible));
end