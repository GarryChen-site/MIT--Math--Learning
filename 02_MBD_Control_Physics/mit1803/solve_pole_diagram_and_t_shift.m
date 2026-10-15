function [pole_diagnostics, sys_models, sim_data] = solve_pole_diagram_and_t_shift(num_vec, den_vec, delay_a, t_span)
% SOLVE_POLE_DIAGRAM_AND_T_SHIFT Analyzes s-domain Pole Diagrams and simulates t-Shift delays.
% Extracts dominant rightmost poles, Region of Convergence (ROC), and long-term asymptotics.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, impulse, pole, zero).
%
% % MBD Concept Mapping: Implements dominant pole analysis (rightmost poles in s-plane)
% %                      and models pure transport delay e^{-as} via LTI InputDelay.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 29 Example 3):
%   num_vec - Numerator polynomial coefficients of F(s) (default: [1, 4, 5] for 1/s + 2/(s^2+2s+5))
%   den_vec - Denominator polynomial coefficients of F(s) (default: [1, 2, 5, 0] for s*(s^2+2s+5))
%   delay_a - Time-shift delay parameter a >= 0 for f_a(t) = u(t-a)f(t-a) (default: 1.5)
%   t_span  - [t_start, t_end] simulation interval (default: [0, 8.0])
%
% Outputs:
%   pole_diagnostics - Struct containing all poles, dominant poles, ROC, and stability regime
%   sys_models       - Struct containing unshifted F(s) and delayed e^{-as}F(s) tf objects
%   sim_data         - Struct containing time vector, unshifted f(t), and delayed f_a(t)

    arguments
        num_vec (1,:) double {mustBeNumeric} = [1.0, 4.0, 5.0]
        den_vec (1,:) double {mustBeNumeric} = [1.0, 2.0, 5.0, 0.0]
        delay_a (1,1) double {mustBeNonnegative} = 1.5
        t_span (1,2) double {mustBeNumeric} = [0, 8.0]
    end

    %% 1. LTI Model Construction with Transport Delay (t-Shift Rule)
    F_unshifted = tf(num_vec, den_vec);
    F_delayed   = tf(num_vec, den_vec, 'InputDelay', delay_a);

    sys_models.F_unshifted = F_unshifted;
    sys_models.F_delayed   = F_delayed;

    %% 2. Pole Diagram & Dominant Rightmost Pole Analysis
    all_poles = pole(F_unshifted);
    all_zeros = zero(F_unshifted);

    if isempty(all_poles)
        roc_sigma = -Inf;
        dom_poles = [];
        stability_str = "Entire Complex Plane (Finite-Duration / Impulse Signal)";
    else
        % Locate rightmost pole real part
        roc_sigma = max(real(all_poles));
        dom_mask  = abs(real(all_poles) - roc_sigma) < 1e-6;
        dom_poles = all_poles(dom_mask);

        if roc_sigma < -1e-9
            stability_str = sprintf('Asymptotically Stable (Exponential Decay ~ e^{%.2ft})', roc_sigma);
        elseif abs(roc_sigma) <= 1e-9
            % Check multiplicity on imaginary axis
            if length(dom_poles) > 1 && abs(dom_poles(1) - dom_poles(2)) < 1e-6
                stability_str = "Unstable Resonance on Imaginary Axis (Polynomial Growth ~ t^k)";
            else
                stability_str = "Marginally Stable / Sub-Exponential (Sustained Steady-State)";
            end
        else
            stability_str = sprintf('Unstable (Exponential Growth ~ e^{+%.2ft})', roc_sigma);
        end
    end

    % Extract oscillatory frequencies from imaginary parts of poles
    osc_freqs = unique(abs(imag(all_poles(abs(imag(all_poles)) > 1e-6))));

    pole_diagnostics.all_poles = all_poles;
    pole_diagnostics.all_zeros = all_zeros;
    pole_diagnostics.dominant_poles = dom_poles;
    pole_diagnostics.roc_abscissa = roc_sigma;
    pole_diagnostics.roc_string = sprintf('Re(s) > %.4f', roc_sigma);
    pole_diagnostics.oscillation_frequencies_rad_s = osc_freqs;
    pole_diagnostics.long_term_behavior = stability_str;
    pole_diagnostics.delay_a = delay_a;

    %% 3. Time-Domain Inversion via Control System Toolbox 'impulse'
    t_eval = linspace(t_span(1), t_span(2), 1000);
    [f_unshifted, ~] = impulse(F_unshifted, t_eval);
    [f_delayed, ~]   = impulse(F_delayed,   t_eval);

    sim_data.t = t_eval;
    sim_data.f_unshifted = f_unshifted;
    sim_data.f_delayed   = f_delayed;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Pole Diagram & t-Shift Solver: MIT 18.03 Lec 29  \n');
    fprintf('====================================================\n');
    fprintf('s-Domain Rational Function F(s):\n');
    display(F_unshifted);
    fprintf('Extracted Poles (Total: %d):\n', length(all_poles));
    for k = 1:length(all_poles)
        fprintf('  p_%d = %+.4f %+.4fi\n', k, real(all_poles(k)), imag(all_poles(k)));
    end
    fprintf('Rightmost Dominant Pole(s) : Re(s) = %+.4f\n', roc_sigma);
    fprintf('Region of Convergence (ROC): %s\n', pole_diagnostics.roc_string);
    fprintf('Long-Term Behavior         : %s\n', stability_str);
    fprintf('t-Shift Parameter (a)      : %.4f s (Factor: e^{-%.2fs}, Poles Unchanged)\n', ...
            delay_a, delay_a);
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [120, 120, 1100, 460]);

        % Left Subplot: s-Plane Pole Diagram & ROC Shading
        subplot(1, 2, 1);
        hold on; grid on; box on;

        x_min = min([real(all_poles); -2.0]) - 1.0;
        x_max = max([real(all_poles);  1.0]) + 1.5;
        y_max = max([abs(imag(all_poles)); 2.0]) + 1.0;

        % Shade ROC half-plane to the right of the rightmost pole
        if ~isinf(roc_sigma)
            patch([roc_sigma, x_max, x_max, roc_sigma], [-y_max, -y_max, y_max, y_max], ...
                  [0.85 0.95 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
                  'DisplayName', sprintf('ROC (%s)', pole_diagnostics.roc_string));
            xline(roc_sigma, 'r--', 'LineWidth', 1.6, ...
                  'DisplayName', sprintf('Rightmost Pole Line Re(s)=%.2f', roc_sigma));
        end

        xline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        yline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');

        if ~isempty(all_poles)
            plot(real(all_poles), imag(all_poles), 'bx', 'MarkerSize', 10, ...
                 'LineWidth', 2.2, 'DisplayName', 'All Poles of F(s)');
            plot(real(dom_poles), imag(dom_poles), 'ro', 'MarkerSize', 10, ...
                 'LineWidth', 2.0, 'DisplayName', 'Dominant Rightmost Pole(s)');
        end
        if ~isempty(all_zeros)
            plot(real(all_zeros), imag(all_zeros), 'go', 'MarkerSize', 8, ...
                 'LineWidth', 1.8, 'DisplayName', 'Zeros of F(s)');
        end

        xlabel('Real Axis Re(s) = \sigma', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Imaginary Axis Im(s) = \omega', 'FontSize', 10, 'FontWeight', 'bold');
        title('s-Plane Pole Diagram & Region of Convergence', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([x_min, x_max]);
        ylim([-y_max, y_max]);

        % Right Subplot: Time-Domain Signals (Unshifted vs. t-Shifted)
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, f_unshifted, 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Unshifted f(t) = L^{-1}\{F(s)\}');
        plot(t_eval, f_delayed, 'r--', 'LineWidth', 2.0, ...
             'DisplayName', sprintf('t-Shifted f_a(t) = L^{-1}\\{e^{-%.1fs}F(s)\\}', delay_a));
        xline(delay_a, 'k:', 'LineWidth', 1.2, ...
              'DisplayName', sprintf('Delay Onset t = %.1fs', delay_a));

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Signal Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('t-Domain: %s', stability_str), 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end