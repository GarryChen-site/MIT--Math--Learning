function report = solve_dynamic_capital_accumulation(borrow_rate_handle, annual_interest_r, total_years_T, num_periods)
% SOLVE_DYNAMIC_CAPITAL_ACCUMULATION Evaluates capital and debt accumulation with interest.
% Continuous formulation:
%   Principal : P = \int_0^T f(t) dt
%   Total Debt: D = \int_0^T f(t) * [ 1 + r*(T - t) ] dt
% Compares continuous integration against discrete periodic repayments.
%
% Syntax:
%   report = solve_dynamic_capital_accumulation(borrow_rate_handle, annual_interest_r, total_years_T, num_periods)
%
% Inputs:
%   borrow_rate_handle - Function handle for borrowing flow rate f(t) in $/year (default = 12000)
%   annual_interest_r  - Annual interest rate r (default = 0.06 for 6%)
%   total_years_T      - Duration of borrowing horizon in years (default = 1.0)
%   num_periods        - Number of discrete periods within T (default = 12 for monthly)
%
% Outputs:
%   report             - Struct containing continuous and discrete debt evaluations

arguments
    borrow_rate_handle  function_handle = @(t) 12000.0 * ones(size(t))
    annual_interest_r   (1,1) double {mustBeNonnegative} = 0.06
    total_years_T       (1,1) double {mustBePositive} = 1.0
    num_periods         (1,1) double {mustBeInteger, mustBePositive} = 12
end

% 1. Continuous Integral Formulation
% Principal integral: \int_0^T f(t) dt
principal_continuous = integral(borrow_rate_handle, 0, total_years_T);

% Total owed integrand: f(t) * (1 + r*(T - t))
total_owed_kernel = @(t) borrow_rate_handle(t) .* (1.0 + annual_interest_r * (total_years_T - t));
total_owed_continuous = integral(total_owed_kernel, 0, total_years_T);
interest_continuous = total_owed_continuous - principal_continuous;

% 2. Discrete Periodic Summation Formulation
dt_period = total_years_T / num_periods;
t_endpoints = (1:num_periods) * dt_period; % Evaluated at the end of each period

discrete_principal_terms = borrow_rate_handle(t_endpoints) * dt_period;
discrete_principal_total = sum(discrete_principal_terms);

% Period interest multiplier: (1 + r * (T - t_i))
discrete_owed_terms = discrete_principal_terms .* (1.0 + annual_interest_r * (total_years_T - t_endpoints));
discrete_owed_total = sum(discrete_owed_terms);
discrete_interest_total = discrete_owed_total - discrete_principal_total;

% Discrepancy between discrete settlement and continuous model
discretization_gap = abs(total_owed_continuous - discrete_owed_total);

% Populate report structure
report.borrowing_duration_years = total_years_T;
report.interest_rate = annual_interest_r;
report.discrete_periods = num_periods;
report.continuous_principal = principal_continuous;
report.continuous_interest = interest_continuous;
report.continuous_total_owed = total_owed_continuous;
report.discrete_principal = discrete_principal_total;
report.discrete_interest = discrete_interest_total;
report.discrete_total_owed = discrete_owed_total;
report.discretization_error = discretization_gap;

% Console report
fprintf('=== Dynamic Capital & Interest Accumulation Analysis ===\n');
fprintf('Horizon Duration (T)       : %.2f Year(s) | Annual Rate r: %.2f%%\n', ...
    total_years_T, annual_interest_r * 100);
fprintf('Discrete Intervals Count   : %d Periods (\\Delta t = %.4f Years)\n', ...
    num_periods, dt_period);
fprintf('Continuous Model:\n');
fprintf('  Net Principal Borrowed   : $%.2f\n', principal_continuous);
fprintf('  Total Interest Accrued   : $%.2f\n', interest_continuous);
fprintf('  Total Repayment Balance  : $%.2f\n', total_owed_continuous);
fprintf('Discrete Summation Model:\n');
fprintf('  Discrete Principal Sum   : $%.2f\n', discrete_principal_total);
fprintf('  Discrete Repayment Total : $%.2f\n', discrete_owed_total);
fprintf('Discretization Residual Gap: $%.4f\n', discretization_gap);
end