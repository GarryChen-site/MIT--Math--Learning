function report = solve_orthogonal_trajectories(slope_expr_sym, test_x, test_y)
% SOLVE_ORTHOGONAL_TRAJECTORIES Derives the orthogonal trajectory differential
% equation and analytical family for a given directional field dy/dx = f(x, y).
% Formula: dy/dx_orthogonal = -1 / f(x, y)
%
% Syntax:
%   report = solve_orthogonal_trajectories(slope_expr_sym, test_x, test_y)
%
% Inputs:
%   slope_expr_sym - Symbolic slope expression for original family (default: 2*y/x for parabolas)
%   test_x         - Arbitrary coordinate x to verify orthogonality (default = 2.0)
%   test_y         - Arbitrary coordinate y to verify orthogonality (default = 2.0)
%
% Outputs:
%   report         - Struct containing original slope, orthogonal ODE, and orthogonal integral

arguments
    slope_expr_sym = []
    test_x         (1,1) double {mustBeNonzero} = 2.0
    test_y         (1,1) double {mustBeNonzero} = 2.0
end

syms x y c_const real
if isempty(slope_expr_sym), slope_expr_sym = 2*y / x; end

% 1. Construct orthogonal slope: m_ortho = -1 / m_orig
orthogonal_slope_sym = simplify(-1.0 / slope_expr_sym);

% 2. Solve orthogonal separable ODE: dy/dx = orthogonal_slope_sym
% Example for -x/(2y): 2y dy = -x dx => y^2 + x^2/2 = c => x^2/4 + y^2/2 = c/2
% Express in separated differential components: M(x)dx + N(y)dy = 0
[num, den] = numden(orthogonal_slope_sym);
% dy/dx = num / den => den * dy - num * dx = 0
integral_curve = simplify(int(den, y) - int(num, x));

% 3. Numerical verification of orthogonality at sample intersection point
m1_val = double(subs(slope_expr_sym, [x, y], [test_x, test_y]));
m2_val = double(subs(orthogonal_slope_sym, [x, y], [test_x, test_y]));
slope_product = m1_val * m2_val;
is_orthogonal = (abs(slope_product - (-1.0)) < 1e-9);

report.original_slope_field = char(slope_expr_sym);
report.orthogonal_slope_field = char(orthogonal_slope_sym);
report.orthogonal_first_integral = [char(integral_curve), ' = c'];
report.test_point = [test_x, test_y];
report.original_slope_at_point = m1_val;
report.orthogonal_slope_at_point = m2_val;
report.slopes_product = slope_product;
report.is_orthogonal_verified = is_orthogonal;

% Console report
fprintf('=== Orthogonal Trajectories Geometric Solver ===\n');
fprintf('Original Curve Slope dy/dx : %s\n', report.original_slope_field);
fprintf('Orthogonal Slope dy/dx     : %s (-1 / m)\n', report.orthogonal_slope_field);
fprintf('Orthogonal First Integral  : %s\n', report.orthogonal_first_integral);
fprintf('Evaluation at Point (%.2f, %.2f):\n', test_x, test_y);
fprintf('  Original Tangent Slope   : %.4f\n', m1_val);
fprintf('  Orthogonal Tangent Slope : %.4f\n', m2_val);
fprintf('  Slope Product (m1 * m2)  : %.4f (Theoretical: -1.0000)\n', slope_product);
fprintf('  Perpendicular Check      : %s\n', mat2str(is_orthogonal));
end