function [fourier_coeffs, sol_analytic, sim_data] = solve_fourier_series_expansion(L, Amp, N_max, sys_plant, t_span)
% SOLVE_FOURIER_SERIES_EXPANSION Calculates Fourier coefficients and simulates LTI filtering.
% Analyzes odd/even parity and reconstructs truncated series: f(t) ~ sum bn*sin(n*pi*t/L).
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, lsim, isproper).
%
% % MBD Concept Mapping: Encapsulates harmonic analysis for non-sinusoidal periodic drives.
% %                      Integrates frequency synthesis with continuous LTI plant G(s)
% %                      to verify actuator filtering and attenuation of high-order harmonics.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 20 Benchmark):
%   L         - Half-period of symmetric window [-L, L] (default: pi)
%   Amp       - Square wave amplitude height (default: pi/4 ~ 0.7854)
%   N_max     - Maximum harmonic truncation index (default: 11)
%   sys_plant - LTI transfer function from Control System Toolbox (default: 2nd-order filter)
%   t_span    - [t_start, t_end] evaluation interval (default: [-pi, 2*pi])
%
% Outputs:
%   fourier_coeffs - Struct containing a0, an, bn vectors and harmonic frequencies
%   sol_analytic   - Struct containing closed-form evaluation function handles
%   sim_data       - Struct containing lsim simulation data [t, u_synth, y_plant]

    arguments
        L (1,1) double {mustBePositive} = pi
        Amp (1,1) double {mustBeNumeric} = pi / 4
        N_max (1,1) double {mustBeInteger, mustBePositive} = 11
        sys_plant = tf(4.0, [1.0, 1.4, 4.0]) % Default 2nd-order plant: wn=2, zeta=0.35
        t_span (1,2) double {mustBeNumeric} = [-pi, 2 * pi]
    end

    %% 1. Analytical Fourier Coefficient Derivation (Odd Square Wave)
    % Parity: f(t) is strictly odd on [-L, L] -> a0 = 0, an = 0 for all n
    w0 = pi / L; % Fundamental circular frequency
    
    n_vec = 1:N_max;
    an_vals = zeros(1, N_max);
    bn_vals = zeros(1, N_max);

    % bn = (2/L) * integral_0^L Amp * sin(n*pi*t/L) dt = (2*Amp/(n*pi)) * [1 - (-1)^n]
    for idx = 1:N_max
        n = n_vec(idx);
        if mod(n, 2) == 1 % Odd harmonic
            bn_vals(idx) = (4 * Amp) / (n * pi);
        else
            bn_vals(idx) = 0.0;
        end
    end

    fourier_coeffs.a0 = 0.0;
    fourier_coeffs.an = an_vals;
    fourier_coeffs.bn = bn_vals;
    fourier_coeffs.harmonics = n_vec;
    fourier_coeffs.omega0 = w0;

    %% 2. Function Handle Reconstruction
    % Truncated synthesis handle
    synth_handle = @(t) evaluate_fourier_sum(t, w0, an_vals, bn_vals, N_max);
    sol_analytic.f_synth = synth_handle;

    %% 3. Numerical Simulation via Control System Toolbox (lsim)
    t_eval = linspace(t_span(1), t_span(2), 1200);
    u_eval = synth_handle(t_eval);

    if isa(sys_plant, 'tf') || isa(sys_plant, 'zpk') || isa(sys_plant, 'ss')
        [y_sim, ~] = lsim(sys_plant, u_eval, t_eval);
    else
        y_sim = u_eval; % Plant bypassed if invalid
    end

    sim_data.t = t_eval;
    sim_data.u = u_eval;
    sim_data.y = y_sim;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Fourier Series Analyzer: MIT 18.03 Lecture 20    \n');
    fprintf('====================================================\n');
    fprintf('Half-Period (L)         : %.4f s -> Fundamental w0 = %.4f rad/s\n', L, w0);
    fprintf('Waveform Amplitude      : %.4f\n', Amp);
    fprintf('Harmonic Truncation (N) : %d\n', N_max);
    fprintf('Parity Assessment       : Odd Function (a0 = 0, an = 0)\n');
    fprintf('First Non-Zero bn Terms :\n');
    for idx = 1:min(5, N_max)
        if bn_vals(idx) > 0
            fprintf('  b_%-2d (Harmonic %2d w0) = %.4f (Exact: %.4f / %d)\n', ...
                    idx, idx, bn_vals(idx), 4*Amp/pi, idx);
        end
    end
    fprintf('Plant Model (G(s))      :\n');
    display(sys_plant);
    fprintf('====================================================\n\n');

    %% 5. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 1050, 480]);

        % Left Subplot: Harmonic Coefficients Spectrum
        subplot(1, 2, 1);
        stem(n_vec, bn_vals, 'filled', 'LineWidth', 1.5, 'MarkerSize', 5);
        grid on; box on;
        xlabel('Harmonic Order n', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Fourier Coefficient b_n', 'FontSize', 10, 'FontWeight', 'bold');
        title('Harmonic Spectrum |b_n| (Odd Harmonic Decay \propto 1/n)', ...
              'FontSize', 11, 'FontWeight', 'bold');
        xlim([0, N_max + 1]);

        % Right Subplot: Time-Domain Synthesis & LTI Output
        subplot(1, 2, 2);
        hold on; grid on; box on;
        
        plot(t_eval, u_eval, 'r--', 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Synthesized Input u(t) (N = %d)', N_max));
        plot(t_eval, y_sim,  'b-',  'LineWidth', 2.0, ...
             'DisplayName', 'Plant Dynamic Output y(t)');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Signal Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title('Time Domain: Actuator Signal Synthesis & Plant Response', ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northwest', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end

%% Local Helper for Fourier Summation
function val = evaluate_fourier_sum(t_in, w0, an, bn, N)
    val = zeros(size(t_in));
    for k = 1:N
        if an(k) ~= 0
            val = val + an(k) * cos(k * w0 * t_in);
        end
        if bn(k) ~= 0
            val = val + bn(k) * sin(k * w0 * t_in);
        end
    end
end