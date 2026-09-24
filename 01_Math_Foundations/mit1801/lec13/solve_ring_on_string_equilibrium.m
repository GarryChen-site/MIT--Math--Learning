function sol = solve_ring_on_string_equilibrium(peg_a, peg_b, string_length, ring_mass, gravity_accel)
% SOLVE_RING_ON_STRING_EQUILIBRIUM Solves for the static mechanical equilibrium
% coordinates (x*, y*) of a sliding ring on a string suspended between (0, 0) and (a, b).
% Minimizes potential energy: y(x) subject to ||p|| + ||p - [a,b]|| = L.
% Theoretical Equations:
%   sin(alpha) = sin(beta) = a / L
%   y* = 0.5 * [ b - sqrt(L^2 - a^2) ]
%   x* = (a / 2) * [ 1 - b / sqrt(L^2 - a^2) ]
%
% Syntax:
%   sol = solve_ring_on_string_equilibrium(peg_a, peg_b, string_length, ring_mass, gravity_accel)
%
% Inputs:
%   peg_a         - Peg 2 horizontal coordinate (a > 0, default = 4.0)
%   peg_b         - Peg 2 vertical coordinate (default = 2.0)
%   string_length - Total string length L (must exceed Euclidean distance, default = 6.0)
%   ring_mass     - Mass of the ring in kg (default = 1.0)
%   gravity_accel - Acceleration due to gravity in m/s^2 (default = 9.81)
%
% Outputs:
%   sol           - Struct containing equilibrium point, angles, and wire tension forces

arguments
    peg_a         (1,1) double {mustBePositive} = 4.0
    peg_b         (1,1) double = 2.0
    string_length (1,1) double {mustBePositive} = 6.0
    ring_mass     (1,1) double {mustBePositive} = 1.0
    gravity_accel (1,1) double {mustBePositive} = 9.81
end

peg_distance = sqrt(peg_a^2 + peg_b^2);
if string_length <= peg_distance
    error('Geometric impossibility: String length L (%.2f) must exceed peg distance (%.2f).', ...
        string_length, peg_distance);
end

if string_length <= abs(peg_a)
    error('Geometric constraint violated: String length L must exceed peg horizontal span a.');
end

% 1. Intermediate root radical
radical_term = sqrt(string_length^2 - peg_a^2);

% 2. Analytic coordinates derived via implicit differentiation
y_star = 0.5 * (peg_b - radical_term);
x_star = (peg_a / 2.0) * (1.0 - (peg_b / radical_term));

% 3. Tension and angle verification
% sin(alpha) = a / L
sin_angle = peg_a / string_length;
angle_alpha_rad = asin(sin_angle);
angle_alpha_deg = rad2deg(angle_alpha_rad);

% Check left and right geometric angles from vertical
dist_left = sqrt(x_star^2 + y_star^2);
dist_right = sqrt((peg_a - x_star)^2 + (peg_b - y_star)^2);

sin_alpha_geo = x_star / dist_left;
sin_beta_geo = (peg_a - x_star) / dist_right;

% Vertical balance of forces: 2 * T * cos(alpha) = m * g => T = mg / (2*cos(alpha))
cos_angle = sqrt(1.0 - sin_angle^2);
tension_force = (ring_mass * gravity_accel) / (2.0 * cos_angle);

% Assemble solution structure
sol.peg1_coords = [0.0, 0.0];
sol.peg2_coords = [peg_a, peg_b];
sol.string_length = string_length;
sol.equilibrium_x = x_star;
sol.equilibrium_y = y_star;
sol.symmetric_angle_deg = angle_alpha_deg;
sol.left_segment_length = dist_left;
sol.right_segment_length = dist_right;
sol.tension_force_newtons = tension_force;
sol.is_angle_equal = (abs(sin_alpha_geo - sin_beta_geo) < 1e-9);

% Console report
fprintf('=== Ring on a String Mechanical Equilibrium Solution ===\n');
fprintf('Peg 1 Location             : (0.00, 0.00)\n');
fprintf('Peg 2 Location             : (%.2f, %.2f)\n', peg_a, peg_b);
fprintf('String Length (L)          : %.4f (Anchor Separation: %.4f)\n', string_length, peg_distance);
fprintf('Equilibrium Coordinates    : (x* = %.4f, y* = %.4f)\n', x_star, y_star);
fprintf('Symmetric Tension Angles   : alpha = beta = %.2f deg (sin = %.4f)\n', ...
    angle_alpha_deg, sin_angle);
fprintf('Equal Angle Check (Proof)  : %s\n', mat2str(sol.is_angle_equal));
fprintf('String Segment Lengths     : L_left = %.4f | L_right = %.4f (Sum = %.4f)\n', ...
    dist_left, dist_right, dist_left + dist_right);
fprintf('Equilibrium String Tension : %.4f N (for mass = %.2f kg)\n', ...
    tension_force, ring_mass);
end