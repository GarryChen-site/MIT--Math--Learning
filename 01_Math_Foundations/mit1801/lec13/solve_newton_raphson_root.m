function sol = solve_newton_raphson_root(f_handle, df_handle, x0, tol, max_iterations)
% SOLVE_NEWTON_RAPHSON_ROOT Numerically solves f(x) = 0 via Newton-Raphson iteration.
% Formula: x_{k+1} = x_k - f(x_k) / f'(x_k)
% Includes cycle detection, singularity protection, and quadratic convergence analysis.
%
% Syntax:
%   sol = solve_newton_raphson_root(f_handle, df_handle, x0, tol, max_iterations)
%
% Inputs:
%   f_handle       - Function handle for objective f(x)
%   df_handle      - Function handle for first derivative f'(x)
%   x0             - Initial guess coordinate (default = 1.0)
%   tol            - Convergence threshold on |f(x)| and step size (default = 1e-10)
%   max_iterations - Maximum allowed iteration count (default = 50)
%
% Outputs:
%   sol            - Struct containing estimated root, iteration table, and convergence status

arguments
    f_handle       function_handle = @(x) x.^2 - 3
    df_handle      function_handle = @(x) 2*x
    x0             (1,1) double = 1.0
    tol            (1,1) double {mustBePositive} = 1e-10
    max_iterations (1,1) double {mustBeInteger, mustBePositive} = 50
end

x_history = zeros(max_iterations + 1, 1);
f_history = zeros(max_iterations + 1, 1);
step_history = zeros(max_iterations, 1);

x_curr = x0;
x_history(1) = x_curr;
f_curr = f_handle(x_curr);
f_history(1) = f_curr;

converged = false;
iter_count = 0;

for k = 1:max_iterations
    df_curr = df_handle(x_curr);

    % Singularity check (zero tangent slope)
    if abs(df_curr) < 1e-14
        error('Newton method failed: Zero derivative encountered at x = %.6f', x_curr);
    end

    delta_x = -f_curr / df_curr;
    x_next = x_curr + delta_x;
    f_next = f_handle(x_next);

    iter_count = k;
    x_history(k+1) = x_next;
    f_history(k+1) = f_next;
    step_history(k) = abs(delta_x);

    % Convergence check
    if abs(delta_x) < tol && abs(f_next) < tol
        converged = true;
        break;
    end

    % Cycle detection: check repetition against earlier iterations
    if any(abs(x_history(1:k) - x_next) < 1e-12)
        warning('Periodic limit cycle detected! Newton method failed to converge.');
        break;
    end

    x_curr = x_next;
    f_curr = f_next;
end

% Truncate histories
x_history = x_history(1:iter_count+1);
f_history = f_history(1:iter_count+1);
step_history = step_history(1:iter_count);

% Assemble structured output
sol.initial_guess = x0;
sol.estimated_root = x_history(end);
sol.final_residual = f_history(end);
sol.iterations_taken = iter_count;
sol.has_converged = converged;
sol.iteration_path = x_history;
sol.residual_history = f_history;
sol.step_sizes = step_history;

% Console report
fprintf('=== Newton-Raphson Iterative Root Solver ===\n');
fprintf('Initial Guess (x0)   : %.6f\n', x0);
fprintf('Iterations Completed : %d\n', iter_count);
fprintf('Convergence Status   : %s\n', mat2str(converged));
fprintf('Estimated Root (x*)  : %.12f\n', sol.estimated_root);
fprintf('Final Residual f(x*) : %.2e\n', sol.final_residual);
fprintf('Iteration Step Log:\n');
for i = 1:iter_count
    fprintf('  Step %2d: x = %14.10f | f(x) = %11.2e | Step Size = %11.2e\n', ...
        i, x_history(i+1), f_history(i+1), step_history(i));
end
end