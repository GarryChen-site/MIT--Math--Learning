function [sol_data, dyn_metrics] = solve_complex_exponential_trajectory(lambda_pole, z0, t_span)
% SOLVE_COMPLEX_EXPONENTIAL_TRAJECTORY Analyzes complex modal trajectories z' = lambda * z.
%
% % Required Toolboxes: 
% %   - Base MATLAB (Core functionality)
%
% % MBD Concept Mapping: Analyzes continuous LTI system pole modal dynamics.
% %                      Maps real part to growth/decay rate and imaginary part to
% %                      angular oscillation frequency.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 6 Question 6.1):
%   lambda_pole - Complex system pole / eigenvalue (default: 1.0 + 2*pi*1i)
%   z0          - Initial complex state z(0) (default: 1.0 + 0i)
%   t_span      - [t_start, t_end] simulation interval (default: [0, 2.0])
%
% Outputs:
%   sol_data    - Struct containing evaluated time, states, and components:
%                 .t, .z, .real_part, .imag_part, .modulus, .phase_rad
%   dyn_metrics - Struct containing frequency, damping rate, and stability topology

    arguments
        lambda_pole (1,1) double = 1.0 + 2 * pi * 1i
        z0 (1,1) double = 1.0 + 0i
        t_span (1,2) double {mustBeNumeric} = [0, 2.0]
    end

    %% 1. Analytical Dynamic Property Extraction
    a_real = real(lambda_pole);
    b_imag = imag(lambda_pole);

    dyn_metrics.growth_rate = a_real;
    dyn_metrics.frequency_rad_s = abs(b_imag);
    dyn_metrics.frequency_hz = abs(b_imag) / (2 * pi);

    % Stability and Phase Plane Classification
    if abs(a_real) < 1e-9
        dyn_metrics.stability_type = "Marginally Stable (Center / Closed Orbit)";
        dyn_metrics.envelope_time_constant = Inf;
    elseif a_real < 0
        dyn_metrics.stability_type = "Asymptotically Stable (Inward Attracting Spiral)";
        dyn_metrics.envelope_time_constant = -1 / a_real; % Decay time constant
    else
        dyn_metrics.stability_type = "Unstable (Outward Expanding Spiral)";
        dyn_metrics.envelope_time_constant = 1 / a_real;  % Growth time constant
    end

    if b_imag > 0
        dyn_metrics.rotation_direction = "Counterclockwise";
    elseif b_imag < 0
        dyn_metrics.rotation_direction = "Clockwise";
    else
        dyn_metrics.rotation_direction = "Pure Real (No Rotation)";
    end

    %% 2. Trajectory Generation via Complex Exponential Closed Form
    t_eval = linspace(t_span(1), t_span(2), 1200);
    z_eval = z0 * exp(lambda_pole * t_eval);

    sol_data.t = t_eval;
    sol_data.z = z_eval;
    sol_data.real_part = real(z_eval);
    sol_data.imag_part = imag(z_eval);
    sol_data.modulus = abs(z_eval);
    sol_data.phase_rad = angle(z_eval);

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Complex Exponential Trajectory Solver (MIT 18.03)\n');
    fprintf('====================================================\n');
    fprintf('Eigenvalue / Pole (\\lambda) : %+.4f %+.4fi\n', a_real, b_imag);
    fprintf('Initial State z(0)       : %+.4f %+.4fi\n', real(z0), imag(z0));
    fprintf('Growth / Decay Rate (a)  : %+.4f\n', a_real);
    fprintf('Oscillation Frequency (b): %.4f rad/s (%.2f Hz)\n', ...
            dyn_metrics.frequency_rad_s, dyn_metrics.frequency_hz);
    fprintf('Rotation Direction       : %s\n', dyn_metrics.rotation_direction);
    fprintf('Phase Topology           : %s\n', dyn_metrics.stability_type);
    fprintf('Envelope Time Constant   : %.4f s\n', dyn_metrics.envelope_time_constant);
    fprintf('====================================================\n\n');

    %% 4. Standalone Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [120, 120, 1050, 480]);

        % Left Plot: Complex Phase Plane
        subplot(1, 2, 1);
        hold on; grid on; box on;
        xline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        plot(sol_data.real_part, sol_data.imag_part, 'b-', 'LineWidth', 1.8, ...
             'DisplayName', 'Complex Trajectory z(t)');
        plot(real(z0), imag(z0), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, ...
             'DisplayName', 'Start Point z(0)');
        plot(sol_data.real_part(end), sol_data.imag_part(end), 'rs', ...
             'MarkerFaceColor', 'r', 'MarkerSize', 7, 'DisplayName', 'End Point z(t_{end})');

        axis equal;
        xlabel('Real Part Re(z)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Imaginary Part Im(z)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Phase Plane: %s (%s)', dyn_metrics.stability_type, dyn_metrics.rotation_direction), ...
              'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);

        % Right Plot: Temporal Projection (Physical Measurable State)
        subplot(1, 2, 2);
        hold on; grid on; box on;
        plot(sol_data.t, sol_data.real_part, 'b-', 'LineWidth', 1.8, ...
             'DisplayName', 'State Response x(t) = Re(z)');
        plot(sol_data.t, sol_data.modulus, 'r--', 'LineWidth', 1.4, ...
             'DisplayName', 'Upper Envelope +|z(t)| = e^{at}');
        plot(sol_data.t, -sol_data.modulus, 'r--', 'LineWidth', 1.4, ...
             'DisplayName', 'Lower Envelope -|z(t)|');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title('Time Domain: Oscillations & Exponential Envelope', 'FontSize', 10.5, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        hold off;
    end
end