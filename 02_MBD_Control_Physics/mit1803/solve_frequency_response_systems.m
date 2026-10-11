function [sys_tf, freq_metrics, sim_data] = solve_frequency_response_systems(sys_type, params, A_in, omega_drive, t_span)
% SOLVE_FREQUENCY_RESPONSE_SYSTEMS Analyzes frequency response and complex gain.
% Supported system types:
%   'dashpot'    : m*x'' + b*x' + k*x = b*y'    (Band-pass, peak gain = 1 at wr = wn)
%   'first_order': x' + k*x = k*y               (Low-pass tide/cooling model, extrema intersection)
%   'suspension' : m*x'' + b*x' + k*x = b*y' + k*y (Quarter-car suspension, gain can exceed 1)
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, bode, lsim, pole, zero).
%
% % MBD Concept Mapping: Encapsulates transfer operators G(s), evaluates Complex Gain H(w) = G(jw),
% %                      and verifies resonant peak invariance and extrema intersection conditions.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 16 Dashpot Benchmark):
%   sys_type    - String identifier: "dashpot", "first_order", or "suspension" (default: "dashpot")
%   params      - Struct with parameters (default: m=1.0, b=1.0, k=4.0)
%   A_in        - Input amplitude (default: 1.0)
%   omega_drive - Input circular driving frequency in rad/s (default: 2.0 rad/s)
%   t_span      - [t_start, t_end] simulation interval (default: [0, 16.0])
%
% Outputs:
%   sys_tf       - Transfer function object G(s) from Control System Toolbox
%   freq_metrics - Struct containing resonant frequency, peak gain, Complex Gain H, and phase lag
%   sim_data     - Struct containing lsim simulation data [t, u, y]

    arguments
        sys_type (1,1) string {mustBeMember(sys_type, ["dashpot", "first_order", "suspension"])} = "dashpot"
        params struct = struct('m', 1.0, 'b', 1.0, 'k', 4.0)
        A_in (1,1) double {mustBeNumeric} = 1.0
        omega_drive (1,1) double {mustBePositive} = 2.0
        t_span (1,2) double {mustBeNumeric} = [0, 16.0]
    end

    %% 1. Construct Transfer Function G(s) via Control System Toolbox
    switch sys_type
        case "dashpot"
            % G(s) = b*s / (m*s^2 + b*s + k)
            m = params.m; b = params.b; k = params.k;
            sys_tf = tf([b, 0], [m, b, k]);
            freq_metrics.omega_r = sqrt(k / m); % Exactly equal to wn
            freq_metrics.peak_gain = 1.0;        % Always 1.0
            
        case "first_order"
            % G(s) = k / (s + k)
            k_val = params.k;
            sys_tf = tf(k_val, [1, k_val]);
            freq_metrics.omega_r = 0.0;          % Low-pass characteristic
            freq_metrics.peak_gain = 1.0;
            
        case "suspension"
            % G(s) = (b*s + k) / (m*s^2 + b*s + k)
            m = params.m; b = params.b; k = params.k;
            sys_tf = tf([b, k], [m, b, k]);
            % Peak frequency for suspension requires numerical peak evaluation
            [mag_peak, ~, w_out] = bode(sys_tf);
            [freq_metrics.peak_gain, max_idx] = max(squeeze(mag_peak));
            freq_metrics.omega_r = w_out(max_idx);
    end

    %% 2. Complex Gain Evaluation at Driving Frequency s = j*omega
    [mag_val, phase_deg] = bode(sys_tf, omega_drive);
    gain_mag = squeeze(mag_val);
    phase_lag_deg = -squeeze(phase_deg); % Convention: phase lag phi = -angle(H)
    phase_lag_rad = phase_lag_deg * (pi / 180);

    % Complex Gain H(w) = G(jw)
    H_complex = gain_mag * exp(-1i * phase_lag_rad);

    freq_metrics.omega_drive = omega_drive;
    freq_metrics.H_complex = H_complex;
    freq_metrics.gain_magnitude = gain_mag;
    freq_metrics.phase_lag_deg = phase_lag_deg;
    freq_metrics.phase_lag_rad = phase_lag_rad;
    freq_metrics.poles = pole(sys_tf);
    freq_metrics.zeros = zero(sys_tf);

    %% 3. Dynamic Time-Domain Simulation via lsim
    t_eval = linspace(t_span(1), t_span(2), 1000);
    u_eval = A_in * cos(omega_drive * t_eval);
    [y_eval, ~] = lsim(sys_tf, u_eval, t_eval);

    sim_data.t = t_eval;
    sim_data.u = u_eval;
    sim_data.y = y_eval;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Frequency Response Analyzer: MIT 18.03 Lec 16    \n');
    fprintf('====================================================\n');
    fprintf('Selected Architecture   : %s\n', sys_type);
    fprintf('Continuous Plant G(s)   : \n');
    display(sys_tf);
    fprintf('Resonant Frequency (wr) : %.4f rad/s\n', freq_metrics.omega_r);
    fprintf('Peak Maximum Gain       : %.4f\n', freq_metrics.peak_gain);
    fprintf('Driving Frequency (w)   : %.4f rad/s\n', omega_drive);
    fprintf('Complex Gain H(w)       : %+.4f %+.4fi\n', real(H_complex), imag(H_complex));
    fprintf('Magnitude Gain |H(w)|   : %.4f\n', gain_mag);
    fprintf('Phase Lag phi = -Arg(H) : %.4f rad (%.2f degrees)\n', ...
            phase_lag_rad, phase_lag_deg);
    fprintf('====================================================\n\n');

    %% 5. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 1050, 480]);

        % Left Subplot: Bode Magnitude Profile
        subplot(1, 2, 1);
        w_plot = logspace(-1, 1.5, 400);
        [mag_curve, ~] = bode(sys_tf, w_plot);
        semilogx(w_plot, squeeze(mag_curve), 'b-', 'LineWidth', 2.0, ...
                 'DisplayName', 'System Gain |G(j\omega)|');
        hold on; grid on; box on;
        
        plot(omega_drive, gain_mag, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, ...
             'DisplayName', sprintf('Operating Point (\\omega = %.2f, Gain = %.3f)', omega_drive, gain_mag));
        if freq_metrics.omega_r > 0
            xline(freq_metrics.omega_r, 'k:', 'LineWidth', 1.2, ...
                  'DisplayName', sprintf('\\omega_r = %.2f rad/s', freq_metrics.omega_r));
        end
        
        xlabel('Frequency \omega (rad/s)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Gain Magnitude |H(\omega)|', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Bode Magnitude: %s Model', sys_type), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);

        % Right Subplot: Dynamic Tracking in Time Domain
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(t_eval, u_eval, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Input Excitation u(t)');
        plot(t_eval, y_eval, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Response Output y(t)');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Signal Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Time Response: \\omega = %.2f rad/s, Lag = %.1f^\\circ', ...
              omega_drive, phase_lag_deg), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'southeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end