function [sol_analytic, sys_tf, sim_data] = solve_fourier_forced_oscillator(omega_n, N_harmonics, t_span)
% SOLVE_FOURIER_FORCED_OSCILLATOR Solves x'' + omega_n^2 * x = sq(t) via Fourier superposition.
% Detects harmonic resonance: wn in {1, 3, 5, ...} excites secular growth,
% whereas wn in {2, 4, 6, ...} remains strictly bounded.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, lsim, pole).
%
% % MBD Concept Mapping: Demonstrates non-sinusoidal periodic drives exciting internal system poles.
% %                      Validates that resonance occurs only when plant poles hit active harmonic lines
% %                      (wn = 1, 3, 5, ... for square wave, but NOT wn = 2, 4, 6, ...).
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 22 Benchmark):
%   omega_n     - Natural frequency of the harmonic oscillator in rad/s (default: 2.0)
%   N_harmonics - Maximum odd harmonic order to include in summation (default: 15)
%   t_span      - [t_start, t_end] evaluation time interval (default: [0, 20.0])
%
% Outputs:
%   sol_analytic - Struct containing harmonic coefficients, resonance flags, and closed-form handles
%   sys_tf       - Transfer function object G(s) from Control System Toolbox
%   sim_data     - Struct containing lsim simulation data [t, u_input, y_sim]

    arguments
        omega_n (1,1) double {mustBePositive} = 2.0
        N_harmonics (1,1) double {mustBeInteger, mustBePositive} = 15
        t_span (1,2) double {mustBeNumeric} = [0, 20.0]
    end

    %% 1. Harmonic Resonance Detection
    % Square wave contains only odd harmonics: n = 1, 3, 5, ...
    is_odd_integer = (abs(omega_n - round(omega_n)) < 1e-6) && (mod(round(omega_n), 2) == 1);
    resonant_n = round(omega_n);

    sol_analytic.omega_n = omega_n;
    sol_analytic.is_resonant = is_odd_integer;
    sol_analytic.resonant_order = [];

    if is_odd_integer
        sol_analytic.resonant_order = resonant_n;
    end

    %% 2. Analytical Superposition Coefficients (Superposition III)
    odd_orders = 1:2:N_harmonics;
    harmonic_weights = zeros(size(odd_orders));

    for idx = 1:length(odd_orders)
        n = odd_orders(idx);
        if is_odd_integer && (n == resonant_n)
            harmonic_weights(idx) = Inf; % Pure resonance singularity
        else
            % Steady-state particular amplitude: B_n = (4 / (pi * n)) / (omega_n^2 - n^2)
            harmonic_weights(idx) = (4 / (pi * n)) / (omega_n^2 - n^2);
        end
    end

    sol_analytic.odd_orders = odd_orders;
    sol_analytic.harmonic_weights = harmonic_weights;

    % Closed-form evaluation handle
    sol_analytic.x_p_handle = @(t) evaluate_forced_response(t, omega_n, odd_orders, harmonic_weights, is_odd_integer, resonant_n);

    %% 3. Control System Toolbox Plant Modeling & Simulation
    % For numerical stability when simulating near-undamped plants, add epsilon damping
    if is_odd_integer
        zeta_sim = 0.002; % Extremely light damping to simulate resonance peak safely
    else
        zeta_sim = 0.0;
    end
    
    sys_tf = tf(1.0, [1.0, 2 * zeta_sim * omega_n, omega_n^2]);

    t_eval = linspace(t_span(1), t_span(2), 1200);
    
    % Synthesize square wave input
    u_eval = zeros(size(t_eval));
    for n = odd_orders
        u_eval = u_eval + (4 / (n * pi)) * sin(n * t_eval);
    end

    [y_sim, ~] = lsim(sys_tf, u_eval, t_eval);

    sim_data.t = t_eval;
    sim_data.u = u_eval;
    sim_data.y = y_sim;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Harmonic Oscillator Solver: MIT 18.03 Lec 22    \n');
    fprintf('====================================================\n');
    fprintf('Natural Frequency (w_n) : %.4f rad/s\n', omega_n);
    fprintf('Harmonic Truncation (N) : %d odd harmonics\n', N_harmonics);
    fprintf('Resonance Condition     : %d\n', sol_analytic.is_resonant);
    
    if sol_analytic.is_resonant
        fprintf('Resonant Harmonic Line  : n = %d (Odd harmonic matches system pole!)\n', resonant_n);
        fprintf('Dynamic Response Regime : Unbounded Secular Growth via RRF\n');
    else
        fprintf('Dynamic Response Regime : Bounded Multi-Harmonic Periodic Orbit\n');
        if abs(omega_n - round(omega_n)) < 1e-6 && mod(round(omega_n), 2) == 0
            fprintf('Special Feature         : w_n is EVEN (%d) -> Input has ZERO energy at this pole!\n', round(omega_n));
        end
    end
    
    fprintf('Transfer Function G(s)  :\n');
    display(sys_tf);
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 1050, 480]);

        % Left Subplot: Harmonic Gain Spectrum
        subplot(1, 2, 1);
        finite_idx = ~isinf(harmonic_weights);
        stem(odd_orders(finite_idx), abs(harmonic_weights(finite_idx)), 'filled', ...
             'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', '|B_n| Harmonic Gain');
        hold on; grid on; box on;
        
        if is_odd_integer
            xline(resonant_n, 'r--', 'LineWidth', 1.5, ...
                  'DisplayName', sprintf('Resonance Line (n = %d)', resonant_n));
        end
        
        xlabel('Harmonic Order n', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Harmonic Weight |B_n|', 'FontSize', 10, 'FontWeight', 'bold');
        title('Harmonic Spectrum Distribution in System Response', 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);

        % Right Subplot: Time-Domain Dynamic Response
        subplot(1, 2, 2);
        hold on; grid on; box on;

        plot(t_eval, u_eval, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0, ...
             'DisplayName', 'Square Input sq(t)');
        plot(t_eval, y_sim, 'b-', 'LineWidth', 1.8, ...
             'DisplayName', sprintf('Oscillator Response x(t), \\omega_n = %.1f', omega_n));

        if is_odd_integer
            secular_scale = (4 / (resonant_n * pi)) / (2 * resonant_n);
            plot(t_eval, secular_scale * t_eval, 'k--', 'LineWidth', 1.2, ...
                 'DisplayName', 'Secular Amplitude Envelope');
            plot(t_eval, -secular_scale * t_eval, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        end

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Displacement x(t)', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Dynamic Tracking: Resonance = %d (\\omega_n = %.1f rad/s)', ...
              sol_analytic.is_resonant, omega_n), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end

%% Local Helper Function for Analytical Superposition
function x_val = evaluate_forced_response(t_in, wn, orders, weights, is_res, res_n)
    x_val = zeros(size(t_in));
    for idx = 1:length(orders)
        n = orders(idx);
        if is_res && (n == res_n)
            % Resonant term: x_res(t) = - (4 / (n*pi)) * (t / (2*n)) * cos(n*t)
            x_val = x_val - (4 / (n * pi)) .* (t_in ./ (2 * n)) .* cos(n * t_in);
        else
            x_val = x_val + weights(idx) .* sin(n * t_in);
        end
    end
end