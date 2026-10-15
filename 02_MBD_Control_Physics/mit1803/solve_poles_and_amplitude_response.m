function [unification_meta, W_tf, freq_data] = solve_poles_and_amplitude_response(p_coeffs, omega_drive, w_max)
% SOLVE_POLES_AND_AMPLITUDE_RESPONSE Unifies s-Plane Pole Diagrams, Transfer Functions & Gain.
% Evaluates W(s) = 1/p(s), extracts the imaginary-axis slice |W(jw)|, and checks stability.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, freqresp, pole, impulse).
%
% % MBD Concept Mapping: Implements the grand unification of MIT 18.03 Lecture 30:
% %                      w(t) --L--> W(s) = 1/p(s) --Slice s=jw--> Complex Gain W(jw).
% %                      Identifies resonant peaks from pole proximity to the imaginary axis.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 30 Underdamped Benchmark):
%   p_coeffs    - Characteristic polynomial coefficients [a_n, ..., a_0] (default: [1.0, 0.8, 4.0])
%   omega_drive - Specific driving circular frequency for sinusoidal steady-state check (default: 2.0)
%   w_max       - Maximum frequency range [-w_max, w_max] for amplitude curve (default: 4.0)
%
% Outputs:
%   unification_meta - Struct with poles, stability flag, complex gain W(jw), and peak locations
%   W_tf             - Transfer function object W(s) = 1/p(s) from Control System Toolbox
%   freq_data        - Struct with two-sided frequency vector, gain curve |W(jw)|, and phase lag

    arguments
        p_coeffs (1,:) double {mustBeNumeric} = [1.0, 0.8, 4.0]
        omega_drive (1,1) double {mustBePositive} = 2.0
        w_max (1,1) double {mustBePositive} = 4.0
    end

    %% 1. Transfer Function W(s) & Stability Analysis via Control System Toolbox
    W_tf = tf(1.0, p_coeffs);
    sys_poles = pole(W_tf);
    rightmost_re = max(real(sys_poles));

    unification_meta.poles = sys_poles;
    unification_meta.rightmost_real_part = rightmost_re;
    unification_meta.is_asymptotically_stable = (rightmost_re < -1e-9);

    if unification_meta.is_asymptotically_stable
        unification_meta.stability_status = "Stable: All poles in Left Half-Plane (All homogeneous solutions are transients)";
    elseif abs(rightmost_re) <= 1e-9
        unification_meta.stability_status = "Marginally Stable / Unstable on Imaginary Axis (Persistent oscillations)";
    else
        unification_meta.stability_status = sprintf('Unstable: Rightmost pole has Re(s) = %+.3f (Grows like e^{%.2ft})', ...
                                                    rightmost_re, rightmost_re);
    end

    %% 2. Imaginary Axis Restriction: Complex Gain W(j*omega)
    w_vec = linspace(-w_max, w_max, 801);
    H_vec = squeeze(freqresp(W_tf, w_vec));
    gain_vec = abs(H_vec);
    phase_lag_vec = -angle(H_vec);

    % Evaluate at specific driving frequency omega_drive
    W_jw_drive = squeeze(freqresp(W_tf, omega_drive));
    unification_meta.omega_drive = omega_drive;
    unification_meta.complex_gain_at_drive = W_jw_drive;
    unification_meta.gain_at_drive = abs(W_jw_drive);
    unification_meta.phase_lag_rad_at_drive = -angle(W_jw_drive);
    unification_meta.phase_lag_deg_at_drive = -angle(W_jw_drive) * (180 / pi);

    % Detect resonant humps on two-sided amplitude response curve
    [pks, locs] = findpeaks(gain_vec, w_vec);
    unification_meta.resonant_peak_gains = pks;
    unification_meta.resonant_peak_freqs = locs;
    unification_meta.num_humps = length(pks);

    freq_data.w_vec = w_vec;
    freq_data.gain_vec = gain_vec;
    freq_data.phase_lag_rad = phase_lag_vec;

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Laplace-Gain Unification: MIT 18.03 Lecture 30   \n');
    fprintf('====================================================\n');
    fprintf('Transfer Function W(s) = 1/p(s):\n');
    display(W_tf);
    fprintf('Stability Assessment    : %s\n', unification_meta.stability_status);
    fprintf('Number of Resonant Humps: %d (Across [-%.1f, +%.1f] rad/s)\n', ...
            unification_meta.num_humps, w_max, w_max);
    if ~isempty(locs)
        fprintf('Peak Frequencies (rad/s): ');
        fprintf('%+.3f  ', locs);
        fprintf('\n');
    end
    fprintf('--- Sinusoidal Response at omega = %.2f rad/s ---\n', omega_drive);
    fprintf('Complex Gain W(j*w)     : %+.4f %+.4fi\n', real(W_jw_drive), imag(W_jw_drive));
    fprintf('Amplitude Gain |W(j*w)| : %.4f\n', unification_meta.gain_at_drive);
    fprintf('Phase Lag -Arg(W(j*w))  : %.4f rad (%.2f deg)\n', ...
            unification_meta.phase_lag_rad_at_drive, unification_meta.phase_lag_deg_at_drive);
    fprintf('====================================================\n\n');

    %% 4. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [100, 100, 1150, 480]);

        % Left Subplot: s-Plane Pole Map & Imaginary Axis Evaluation Line
        subplot(1, 2, 1);
        hold on; grid on; box on;
        xline(0, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Imaginary Axis s = j\omega (Gain Slice)');
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        plot(real(sys_poles), imag(sys_poles), 'bx', 'MarkerSize', 11, ...
             'LineWidth', 2.2, 'DisplayName', 'Poles of W(s)');
        plot([0, 0], [omega_drive, -omega_drive], 'ro', 'MarkerFaceColor', 'r', ...
             'MarkerSize', 7, 'DisplayName', sprintf('Input Poles \\pm j%.1f', omega_drive));

        % Draw distance vectors from poles to j*omega_drive
        for k = 1:length(sys_poles)
            plot([real(sys_poles(k)), 0], [imag(sys_poles(k)), omega_drive], ...
                 'k--', 'LineWidth', 0.9, 'HandleVisibility', 'off');
        end

        xlabel('Real Axis \sigma = Re(s)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Imaginary Axis \omega = Im(s)', 'FontSize', 10, 'FontWeight', 'bold');
        title('s-Plane Pole Diagram & Imaginary Axis Proximity', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([min([real(sys_poles); -2.0]) - 1.0, max([real(sys_poles); 0.5]) + 1.0]);
        ylim([-w_max, w_max]);

        % Right Subplot: Two-Sided Amplitude Response Curve |W(j*omega)|
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(w_vec, gain_vec, 'b-', 'LineWidth', 2.0, ...
             'DisplayName', 'Amplitude Response |W(j\omega)|');
        if ~isempty(locs)
            plot(locs, pks, 'mv', 'MarkerFaceColor', 'm', 'MarkerSize', 7, ...
                 'DisplayName', sprintf('Resonant Humps (%d peaks)', length(locs)));
        end
        plot(omega_drive, unification_meta.gain_at_drive, 'ro', 'MarkerFaceColor', 'r', ...
             'MarkerSize', 7, 'DisplayName', sprintf('Gain at \\omega=%.1f: %.3f', ...
             omega_drive, unification_meta.gain_at_drive));

        xline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        xlabel('Circular Frequency \omega (rad/s)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Gain Magnitude |W(j\omega)|', 'FontSize', 10, 'FontWeight', 'bold');
        title('Imaginary-Axis Slice: Two-Sided Amplitude Response', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([-w_max, w_max]);
        hold off;
    end
end