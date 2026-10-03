function sol = solve_block_stacking_equilibrium(block_length, block_thickness, num_blocks, target_overhang)
% SOLVE_BLOCK_STACKING_EQUILIBRIUM Computes the maximum stable overhang
% for a stack of identical rigid blocks using recursive center-of-mass balancing (Lecture 36).
%
% Formulas:
%   Recursive Offset : delta_x_k = (block_length / 2) * (1 / k),  for k = 1 to N-1
%   Total Overhang   : S_overhang = (block_length / 2) * \sum_{k=1}^{N-1} (1 / k)
%   Inverse Estimate : Required N approx exp( 2 * D_target / block_length )
%
% Syntax:
%   sol = solve_block_stacking_equilibrium(block_length, block_thickness, num_blocks, target_overhang)
%
% Inputs:
%   block_length    - Length of single block L in meters (default = 0.30 m)
%   block_thickness - Thickness/height of single block H in meters (default = 0.03 m)
%   num_blocks      - Total number of stacked blocks N (default = 6)
%   target_overhang - Optional target horizontal overhang distance to reach (default = 0.60 m)
%
% Outputs:
%   sol             - Struct containing block coordinates, total overhang, and inverse estimation

arguments
    block_length    (1,1) double {mustBePositive} = 0.30
    block_thickness (1,1) double {mustBePositive} = 0.03
    num_blocks      (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(num_blocks, 2)} = 6
    target_overhang (1,1) double {mustBePositive} = 0.60
end

half_L = block_length / 2.0;

% 1. Stepwise Overhang Offsets from Top Down: delta_k = (L/2) / k
k_vec = 1:(num_blocks - 1);
delta_offsets = half_L ./ k_vec;
harmonic_sum = sum(1.0 ./ k_vec);
total_overhang_current = half_L * harmonic_sum;

% 2. Coordinate Reconstruction (Bottom block at x = 0, y = 0)
x_left_coords = zeros(num_blocks, 1);
y_coords = zeros(num_blocks, 1);

for i = (num_blocks - 1):-1:1
    x_left_coords(i) = x_left_coords(i + 1) + delta_offsets(i);
    y_coords(i) = (num_blocks - i) * block_thickness;
end

x_right_coords = x_left_coords + block_length;

% 3. Verify Recursive Center of Mass Stability
% CoM of top m blocks must not exceed right edge of (m+1)-th block
is_stable = true;
com_margins = zeros(num_blocks - 1, 1);
for m = 1:(num_blocks - 1)
    % Mass is uniform => average of individual block centers
    com_top_m = mean(x_left_coords(1:m) + half_L);
    edge_supporting = x_right_coords(m + 1);
    margin = edge_supporting - com_top_m;
    com_margins(m) = margin;
    if margin < -1e-9
        is_stable = false;
    end
end

% 4. Inverse Problem: Blocks required to span target_overhang D
% D = (L/2) * \sum_{k=1}^N (1/k) approx (L/2) * ln(N) => N approx exp(2*D / L)
required_harmonic_val = target_overhang / half_L;
% Exact search via ceiling of exponential approximation
N_estimated_rough = exp(required_harmonic_val);

% Check realistic engineering height
stack_height_estimated_meters = N_estimated_rough * block_thickness;

% Assemble structured output
sol.block_length_m = block_length;
sol.block_thickness_m = block_thickness;
sol.num_blocks = num_blocks;
sol.overhang_offsets = delta_offsets';
sol.total_overhang_meters = total_overhang_current;
sol.overhang_in_block_lengths = total_overhang_current / block_length;
sol.is_static_equilibrium_held = is_stable;
sol.min_stability_margin = min(com_margins);
sol.target_overhang_m = target_overhang;
sol.blocks_required_for_target = N_estimated_rough;
sol.height_required_for_target_m = stack_height_estimated_meters;

% Console report
fprintf('=== Multi-Body Block Stacking Harmonic Equilibrium ===\n');
fprintf('Block Specs                 : Length = %.3f m | Thickness = %.3f m\n', block_length, block_thickness);
fprintf('Stack Configuration         : %d Blocks (Total Height = %.3f m)\n', num_blocks, num_blocks * block_thickness);
fprintf('Current Achieved Overhang   : %.4f m (%.2f%% of single block length)\n', ...
    total_overhang_current, sol.overhang_in_block_lengths * 100);
fprintf('Static Stability Status     : Confirmed Equilibrium: %s (Zero Tip Margin: %.2e m)\n', ...
    mat2str(is_stable), sol.min_stability_margin);
fprintf('--- Inverse Reach Solver ---\n');
fprintf('Target Horizontal Overhang  : %.4f m (%.2f Block Lengths)\n', target_overhang, target_overhang / block_length);
fprintf('Required Block Count (N)    : ~%.2e Blocks (Theoretical: e^{%.2f})\n', ...
    N_estimated_rough, required_harmonic_val);
fprintf('Resulting Stack Height      : ~%.2e meters\n', stack_height_estimated_meters);
end