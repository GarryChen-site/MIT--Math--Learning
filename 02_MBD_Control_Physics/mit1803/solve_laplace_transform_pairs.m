function [F_tf, laplace_meta, sim_data] = solve_laplace_transform_pairs(poly_coeffs, a_exp, A_cos, B_sin, omega, t_span)
% SOLVE_LAPLACE_TRANSFORM_PAIRS Synthesizes s-domain rational function F(s) = L{f(t)}.
% Signal form: f(t) = exp(a_exp*t) * [ sum_{k=0}^N c_k * t^k + A_cos*cos(omega*t) + B_sin*sin(omega*t) ]
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, impulse, pole, zero, bode).
%
% % MBD Concept Mapping: Implements linearity, power rule L[t^n] = n!/s^{n+1}, and s-shift F(s - a).
% %                      Constructs rational LTI transfer function F(s), identifies Region of Convergence (ROC),
% %                      and validates time-domain reconstruction via impulse response kernel equivalence.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 26 Question 26.1: (1+t)^2):
%   poly_coeffs - Vector [c0, c1, ..., cN] for c0 + c1*t + ... + cN*t^N (default: [1, 2, 1] for (1+t)^2)
%   a_exp       - Real exponential shift parameter in exp(a_exp*t) (default: 0.0)
%   A_cos       - Cosine component amplitude (default: 0.0)
%   B_sin       - Sine component amplitude (default: 0.0)
%   omega       - Sinusoidal circular frequency in rad/s (default: 1.0)
%   t_span      - [t_start, t_end] time evaluation interval (default: [0, 4.0])
%
% Outputs:
%   F_tf         - Rational transfer function object F(s) from Control System Toolbox
%   laplace_meta - Struct containing ROC boundary, poles, zeros, and verification error
%   sim_data     - Struct containing time vector, analytical f(t), and impulse-recovered f_inv(t)

    arguments
        poly_coeffs (1,:) double {mustBeNumeric} = [1.0, 2.0, 1.0] % (1 + t)^2 = 1 + 2*t + 1*t^2
        a_exp (1,1) double {mustBeNumeric} = 0.0
        A_cos (1,1) double {mustBeNumeric} = 0.0
        B_sin (1,1) double {mustBeNumeric} = 0.0
        omega (1,1) double {mustBePositive} = 1.0
        t_span (1,2) double {mustBeNumeric} = [0, 4.0]
    end

    %% 1. Symbolic-Free Rational Synthesis via Control System Toolbox 'tf'
    s = tf('s');
    s_shifted = s - a_exp; % Apply Rule 2 (s-shift): s -> (s - a)

    F_total = tf(0, 1);

    % Part A: Polynomial powers sum_{k=0}^N c_k * (k!) / (s - a)^{k+1}
    for idx = 1:length(poly_coeffs)
        n_pow = idx - 1;
        c_val = poly_coeffs(idx);
        if abs(c_val) > 0
            F_term = (c_val * factorial(n_pow)) / (s_shifted^(n_pow + 1));
            F_total = F_total + F_term;
        end
    end

    % Part B: Sinusoidal components with s-shift
    if abs(A_cos) > 0
        F_cos = A_cos * (s_shifted / (s_shifted^2 + omega^2));
        F_total = F_total + F_cos;
    end
    if abs(B_sin) > 0
        F_sin = B_sin * (omega / (s_shifted^2 + omega^2));
        F_total = F_total + F_sin;
    end

    % Simplify rational transfer function representation
    F_tf = minreal(F_total, 1e-9);

    %% 2. Region of Convergence (ROC) & Pole-Zero Extraction
    sys_poles = pole(F_tf);
    sys_zeros = zero(F_tf);
    if isempty(sys_poles)
        roc_bound = -Inf;
    else
        roc_bound = max(real(sys_poles));
    end

    laplace_meta.poles = sys_poles;
    laplace_meta.zeros = sys_zeros;
    laplace_meta.roc_abscissa = roc_bound;
    laplace_meta.roc_string = sprintf('Re(s) > %.4f', roc_bound);

    %% 3. Time-Domain Verification: Analytical vs. Control Toolbox 'impulse'
    t_eval = linspace(t_span(1), t_span(2), 800);
    
    % Exact analytical evaluation of f(t)
    poly_val = zeros(size(t_eval));
    for idx = 1:length(poly_coeffs)
        n_pow = idx - 1;
        poly_val = poly_val + poly_coeffs(idx) * (t_eval.^n_pow);
    end
    trig_val = A_cos * cos(omega * t_eval) + B_sin * sin(omega * t_eval);
    f_analytic = exp(a_exp * t_eval) .* (poly_val + trig_val);

    % Inverse Laplace via Control System Toolbox impulse response
    [f_impulse_inv, ~] = impulse(F_tf, t_eval);

    max_err = max(abs(f_analytic(:) - f_impulse_inv(:)));
    laplace_meta.verification_max_error = max_err;

    sim_data.t = t_eval;
    sim_data.f_analytic = f_analytic;
    sim_data.f_recovered = f_impulse_inv;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Laplace Transform Solver: MIT 18.03 Lecture 26   \n');
    fprintf('====================================================\n');
    fprintf('s-Domain Rational Function F(s) :\n');
    display(F_tf);
    fprintf('Region of Convergence (ROC)     : %s\n', laplace_meta.roc_string);
    fprintf('Inverse LT Verification Error   : %e (|f(t) - impulse(F(s))|)\n', max_err);
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 1050, 460]);

        % Left Subplot: Time-Domain Signal Reconstruction
        subplot(1, 2, 1);
        hold on; grid on; box on;
        plot(t_eval, f_analytic, 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Analytical Signal f(t) = (1 + t)^2');
        plot(t_eval(1:30:end), f_impulse_inv(1:30:end), 'ro', 'MarkerSize', 6, ...
             'DisplayName', 'Recovered L^{-1}\{F(s)\} via impulse(F)');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Signal Amplitude f(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title('t-Domain Verification: f(t) vs. L^{-1}\{F(s)\}', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);

        % Right Subplot: Pole-Zero Map & ROC in the s-Plane
        subplot(1, 2, 2);
        hold on; grid on; box on;
        xline(roc_bound, 'r--', 'LineWidth', 1.5, ...
              'DisplayName', sprintf('ROC Boundary: Re(s) = %.2f', roc_bound));
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        if ~isempty(sys_poles)
            plot(real(sys_poles), imag(sys_poles), 'rx', 'MarkerSize', 10, ...
                 'LineWidth', 2.0, 'DisplayName', sprintf('Poles (%d total)', length(sys_poles)));
        end
        if ~isempty(sys_zeros)
            plot(real(sys_zeros), imag(sys_zeros), 'bo', 'MarkerSize', 8, ...
                 'LineWidth', 1.8, 'DisplayName', sprintf('Zeros (%d total)', length(sys_zeros)));
        end

        xlabel('Real Axis Re(s)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Imaginary Axis Im(s)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('s-Domain Pole-Zero Map (%s)', laplace_meta.roc_string), ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([min([real(sys_poles); -2]) - 1, max([real(sys_poles); 2]) + 1]);
        ylim([-3, 3]);
        hold off;
    end
end