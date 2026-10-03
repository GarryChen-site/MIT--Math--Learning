function report = solve_matrix_exponential_taylor(A_mat, truncation_order_N, time_horizon_tf)
% SOLVE_MATRIX_EXPONENTIAL_TAYLOR Computes the state transition matrix
% Phi(t) = e^{At} for continuous LTI systems via matrix Taylor expansion:
% Formula: \Phi_N(t) = \sum_{k=0}^N \frac{(A \cdot t)^k}{k!} (Lecture 37 Example 7 generalization).
% Analyzes matrix-norm convergence against benchmark MATLAB expm(A*t).
%
% Syntax:
%   report = solve_matrix_exponential_taylor(A_mat, truncation_order_N, time_horizon_tf)
%
% Inputs:
%   A_mat              - Square continuous state dynamics matrix (n x n) (default: [0, 1; -4, 0])
%   truncation_order_N - Taylor expansion truncation order (default = 8)
%   time_horizon_tf    - Target evaluation time horizon in seconds (default = 1.0)
%
% Outputs:
%   report             - Struct containing exact/approximate matrices, Frobenius error, and spectral norm

arguments
    A_mat              (:,:) double = [0.0, 1.0; -4.0, 0.0]
    truncation_order_N (1,1) double {mustBeInteger, mustBePositive} = 8
    time_horizon_tf    (1,1) double {mustBePositive} = 1.0
end

[rows, cols] = size(A_mat);
if rows ~= cols
    error('System matrix A must be square.');
end

% 1. Benchmark State Transition via expm
Phi_exact = expm(A_mat * time_horizon_tf);

% 2. Matrix Taylor Series Computation: \sum_{k=0}^N (A*tf)^k / k!
At = A_mat * time_horizon_tf;
Phi_taylor = eye(rows);
At_power = eye(rows);

for k = 1:truncation_order_N
    At_power = At_power * At;
    Phi_taylor = Phi_taylor + At_power / factorial(k);
end

% 3. Norm Discrepancy & Truncation Bounds
matrix_residual = Phi_exact - Phi_taylor;
frobenius_error = norm(matrix_residual, 'fro');
spectral_error = norm(matrix_residual, 2);

% Theoretical Cauchy-Hadamard remainder upper bound:
% || R_N || <= (||At||^(N+1) / (N+1)!) * exp(||At||)
norm_At = norm(At, 2);
tail_bound_theoretical = ((norm_At^(truncation_order_N + 1)) / factorial(truncation_order_N + 1)) * exp(norm_At);

% Assemble structured report
report.system_matrix_A = A_mat;
report.dimension = rows;
report.truncation_order = truncation_order_N;
report.evaluation_time_tf = time_horizon_tf;
report.state_transition_exact = Phi_exact;
report.state_transition_taylor = Phi_taylor;
report.frobenius_norm_error = frobenius_error;
report.spectral_norm_error = spectral_error;
report.theoretical_upper_bound = tail_bound_theoretical;

% Console report
fprintf('=== State Transition Matrix Exponential via Taylor Series ===\n');
fprintf('Matrix Dimension (n x n)   : %d x %d | Horizon tf = %.2f s\n', rows, cols, time_horizon_tf);
fprintf('Taylor Series Truncation   : Order N = %d\n', truncation_order_N);
fprintf('Exact expm(A*tf) Matrix    :\n'); disp(Phi_exact);
fprintf('Taylor Series \\Phi_N(tf)   :\n'); disp(Phi_taylor);
fprintf('Spectral 2-Norm Error      : %.2e (Frobenius: %.2e)\n', spectral_error, frobenius_error);
fprintf('Theoretical Cauchy Bound   : %.2e\n', tail_bound_theoretical);
end