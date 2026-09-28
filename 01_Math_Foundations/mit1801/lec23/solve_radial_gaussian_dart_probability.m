function report = solve_radial_gaussian_dart_probability(p_target, inner_mult, outer_mult, sector_fraction)
% SOLVE_RADIAL_GAUSSIAN_DART_PROBABILITY Solves 2D Gaussian radial probability:
% PDF Kernel       : p(r) = 2 * r * exp(-r^2)  (after angular 2*pi normalization)
% Calibration      : P(0 <= r <= r1) = p_target => 1 - exp(-r1^2) = p_target => exp(-r1^2) = 1 - p_target
% Sector Probability: P_sector = sector_fraction * \int_{inner}^{outer} 2*r*exp(-r^2) dr
% Gaussian Factor  : Validates \int_{-inf}^{inf} exp(-x^2) dx = sqrt(pi).
%
% Syntax:
%   report = solve_radial_gaussian_dart_probability(p_target, inner_mult, outer_mult, sector_fraction)
%
% Inputs:
%   p_target        - Probability of hitting inner dartboard 0 <= r <= r1 (default = 0.5)
%   inner_mult      - Multiplier for sector inner radius r_in = inner_mult * r1 (default = 2.0)
%   outer_mult      - Multiplier for sector outer radius r_out = outer_mult * r1 (default = 3.0)
%   sector_fraction - Circumferential angle fraction (default = 1/6 for 2 hours on a clock)
%
% Outputs:
%   report          - Struct containing calibrated r1, sector probability, and Gaussian proof

arguments
    p_target        (1,1) double {mustBeInRange(p_target, 0, 1)} = 0.5
    inner_mult      (1,1) double {mustBeNonnegative} = 2.0
    outer_mult      (1,1) double {mustBePositive} = 3.0
    sector_fraction (1,1) double {mustBeInRange(sector_fraction, 0, 1)} = 1.0 / 6.0
end

if inner_mult >= outer_mult
    error('Inner radius multiplier must be strictly less than outer multiplier.');
end

% 1. Radial Calibration: P(0 <= r <= r1) = \int_0^r1 2*r*exp(-r^2) dr = 1 - exp(-r1^2) = p_target
% exp(-r1^2) = 1 - p_target => r1^2 = -ln(1 - p_target)
decay_base = 1.0 - p_target; % For p_target = 0.5, decay_base = 0.5
r1_calibrated = sqrt(-log(decay_base));

% 2. Sector Boundaries
r_inner = inner_mult * r1_calibrated;
r_outer = outer_mult * r1_calibrated;

% 3. Radial Ring Probability: \int_{r_in}^{r_out} 2*r*exp(-r^2) dr = exp(-r_in^2) - exp(-r_out^2)
% Note: exp(-r_in^2) = exp(-(inner_mult*r1)^2) = (exp(-r1^2))^(inner_mult^2) = decay_base^(inner_mult^2)
p_ring_analytic = (decay_base^(inner_mult^2)) - (decay_base^(outer_mult^2));

% Multiply by angular sector fraction (e.g., 1/6)
p_sector_analytic = sector_fraction * p_ring_analytic;

% 4. Numerical Quadrature Verification
radial_pdf = @(r) 2.0 * r .* exp(-r.^2);
p_ring_num = integral(radial_pdf, r_inner, r_outer, 'RelTol', 1e-12, 'AbsTol', 1e-14);
p_sector_num = sector_fraction * p_ring_num;

% 5. Gaussian Normalization Check: \int_{-inf}^{inf} exp(-x^2) dx == sqrt(pi)
Q_analytic = sqrt(pi);
Q_num = integral(@(x) exp(-x.^2), -Inf, Inf, 'RelTol', 1e-12, 'AbsTol', 1e-14);

% Populate output structure
report.target_calibration_prob = p_target;
report.calibrated_r1 = r1_calibrated;
report.decay_base_factor = decay_base;
report.sector_inner_radius = r_inner;
report.sector_outer_radius = r_outer;
report.sector_angle_fraction = sector_fraction;
report.analytic_sector_probability = p_sector_analytic;
report.numerical_sector_probability = p_sector_num;
report.sector_probability_percent = p_sector_analytic * 100.0;
report.gaussian_integral_analytic = Q_analytic;
report.gaussian_integral_numerical = Q_num;
report.is_gaussian_verified = (abs(Q_analytic - Q_num) < 1e-10);

% Console report
fprintf('=== 2D Radial Gaussian Dartboard & Hazard Sector Solver ===\n');
fprintf('Calibration Condition       : P(r <= r1) = %.2f%% => r1 = %.4f\n', ...
    p_target * 100, r1_calibrated);
fprintf('Hazard Ring Radial Span     : [%.4f, %.4f] (Multipliers: [%.1f, %.1f]*r1)\n', ...
    r_inner, r_outer, inner_mult, outer_mult);
fprintf('Angular Sector Fraction     : %.4f (i.e., 1/%.1f of full ring)\n', ...
    sector_fraction, 1.0 / sector_fraction);
fprintf('Full Annular Ring Hit Prob  : %.6f (%.2f%%)\n', p_ring_analytic, p_ring_analytic * 100);
fprintf('Hazard Sector Strike Prob   : %.6f (%.2f%%)\n', p_sector_analytic, report.sector_probability_percent);
fprintf('Quadrature Error Check      : %.2e\n', abs(p_sector_analytic - p_sector_num));
fprintf('Gaussian Proof: \\int exp(-x^2)dx = %.10f (Theory \\surd\\pi = %.10f) -> Verified: %s\n', ...
    Q_num, Q_analytic, mat2str(report.is_gaussian_verified));
end