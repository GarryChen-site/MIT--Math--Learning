function report = convert_polar_conic_to_cartesian(semi_latus_p, eccentricity_e, trig_alignment)
% CONVERT_POLAR_CONIC_TO_CARTESIAN Converts a general polar conic section
% centered at a focus: r = p / (1 + e*trig(theta)) into Cartesian standard form (Lecture 32 Example 7).
%
% Algebraic Derivation:
%   r + e*r*sin(theta) = p => sqrt(x^2 + y^2) = p - e*y
%   Squaring: x^2 + y^2 = p^2 - 2*p*e*y + e^2*y^2 => x^2 + (1 - e^2)*y^2 + 2*p*e*y = p^2
%
% Syntax:
%   report = convert_polar_conic_to_cartesian(semi_latus_p, eccentricity_e, trig_alignment)
%
% Inputs:
%   semi_latus_p    - Semi-latus rectum parameter p > 0 (default = 1.0)
%   eccentricity_e  - Orbital eccentricity e >= 0 (default = 0.5 for ellipse)
%   trig_alignment  - Direction char: 'sin' (aligned on y-axis) or 'cos' (aligned on x-axis) (default: 'sin')
%
% Outputs:
%   report          - Struct containing conic classification, Cartesian equation string, and parameters

arguments
    semi_latus_p   (1,1) double {mustBePositive} = 1.0
    eccentricity_e (1,1) double {mustBeNonnegative} = 0.5
    trig_alignment (1,:) char {mustBeMember(trig_alignment, {'sin', 'cos'})} = 'sin'
end

p = semi_latus_p;
e = eccentricity_e;

% 1. Conic Type Classification based on Eccentricity
if abs(e) < 1e-12
    conic_type = 'Circle (e = 0)';
elseif e < 1.0
    conic_type = 'Ellipse (0 < e < 1)';
elseif abs(e - 1.0) < 1e-12
    conic_type = 'Parabola (e = 1)';
else
    conic_type = 'Hyperbola (e > 1)';
end

% 2. Derive Standard Cartesian Form Coefficients:
% For alignment 'sin': x^2 + (1 - e^2)*y^2 + 2*p*e*y = p^2
% For alignment 'cos': (1 - e^2)*x^2 + y^2 + 2*p*e*x = p^2
quad_factor = 1.0 - e^2;
linear_factor = 2.0 * p * e;
const_rhs = p^2;

% 3. Extract Geometric Parameters for Ellipse (e < 1)
if e < 1.0
    % Complete the square: quad_factor * (y + shift)^2
    % shift = (linear_factor) / (2 * quad_factor) = (p * e) / (1 - e^2)
    shift_val = (p * e) / quad_factor;
    RHS_completed = const_rhs + quad_factor * (shift_val^2);

    % Semi-major axis a_axis, Semi-minor axis b_axis
    % (y + shift)^2 / a^2 + x^2 / b^2 = 1 (if aligned on sin)
    a_semi_major = sqrt(RHS_completed / quad_factor);
    b_semi_minor = sqrt(RHS_completed);
    focal_c = sqrt(abs(a_semi_major^2 - b_semi_minor^2));
    enclosed_area = pi * a_semi_major * b_semi_minor;
else
    a_semi_major = NaN;
    b_semi_minor = NaN;
    focal_c = NaN;
    enclosed_area = Inf;
end

% Construct Cartesian Equation String
if strcmp(trig_alignment, 'sin')
    cart_str = sprintf('x^2 + %.4f*y^2 %+.4f*y = %.4f', quad_factor, linear_factor, const_rhs);
else
    cart_str = sprintf('%.4f*x^2 + y^2 %+.4f*x = %.4f', quad_factor, linear_factor, const_rhs);
end

% Assemble structured report
report.semi_latus_rectum_p = p;
report.eccentricity_e = e;
report.orientation = trig_alignment;
report.conic_classification = conic_type;
report.cartesian_equation = cart_str;
report.semi_major_axis = a_semi_major;
report.semi_minor_axis = b_semi_minor;
report.focal_distance = focal_c;
report.enclosed_orbital_area = enclosed_area;

% Console report
fprintf('=== Polar Conic Section to Cartesian Orbit Solver ===\n');
fprintf('Polar Model Formulation    : r = %.2f / [ 1 + %.2f*%s(\\theta) ]\n', p, e, trig_alignment);
fprintf('Geometric Classification   : %s\n', conic_type);
fprintf('Cartesian Algebraic Form   : %s\n', cart_str);
if e < 1.0
    fprintf('Semi-Major Axis (a)        : %.4f\n', a_semi_major);
    fprintf('Semi-Minor Axis (b)        : %.4f\n', b_semi_minor);
    fprintf('Focal Separation (c = ae)  : %.4f\n', focal_c);
    fprintf('Enclosed Orbit Area (\\pi ab): %.10f\n', enclosed_area);
end
end