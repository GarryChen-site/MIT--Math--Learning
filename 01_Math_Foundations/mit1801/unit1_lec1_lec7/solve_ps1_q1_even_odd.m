function sol = solve_ps1_q1_even_odd(f_sym)
% SOLVE_PS1_Q1_EVEN_ODD Decomposes an arbitrary scalar function f(x)
% into its symmetric even component f_e(x) and anti-symmetric odd component f_o(x).
% Formulas: f_e(x) = [f(x) + f(-x)] / 2, f_o(x) = [f(x) - f(-x)] / 2
%
% Syntax:
%   sol = solve_ps1_q1_even_odd(f_sym)
%
% Inputs:
%   f_sym - Symbolic expression in variable 'x' (default: (x-1)/(x+1))
%
% Outputs:
%   sol   - Struct containing even component, odd component, and simplified expressions

arguments
    f_sym = []
end

syms x
if isempty(f_sym), f_sym = (x - 1) / (x + 1); end

% Evaluate f(-x)
f_neg_x = subs(f_sym, x, -x);

% Decompose into even and odd components
even_part = simplify((f_sym + f_neg_x) / 2);
odd_part  = simplify((f_sym - f_neg_x) / 2);

% Verification check: f_e(x) + f_o(x) == f(x)
reconstructed_f = simplify(even_part + odd_part);
is_valid = isequal(reconstructed_f, simplify(f_sym));

sol.original_function = char(f_sym);
sol.even_part = char(even_part);
sol.odd_part = char(odd_part);
sol.is_identity_verified = is_valid;

% Console report
fprintf('=== Problem 1: Even and Odd Signal Decomposition ===\n');
fprintf('Input Function f(x) : %s\n', sol.original_function);
fprintf('Even Component f_e  : %s\n', sol.even_part);
fprintf('Odd Component f_o   : %s\n', sol.odd_part);
fprintf('Reconstruction Check: %s\n', mat2str(is_valid));
end