function [fourier_meta, sol_analytic, sim_data] = solve_fourier_shifted_arbitrary_period(L, t_shift, Amp, N_terms, sys_plant, t_span)
% SOLVE_FOURIER_SHIFTED_ARBITRARY_PERIOD Calculates Fourier series for shifted/scaled square waves.
% Evaluates Dirichlet convergence at jump discontinuities and drives an LTI plant via lsim.
%
% % Required Toolboxes: 
% %   - Control System Toolbox (tf, lsim).
%
% % MBD Concept Mapping: Demonstrates non-sinusoidal periodic drives, transport delays,
% %                      and LTI low-pass plant attenuation of high-order harmonics.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 21 Scaling Benchmark):
%   L         - Half-period of the symmetric window [-L, L] (default: 2.0 -> Period = 4.0)
%   t_shift   - Time-shift offset parameter (default: 0.0)
%   Amp       - Square wave amplitude height (default: 1.0)
%   N_terms   - Number of odd harmonic terms to retain (default: 15)
%   sys_plant - LTI plant transfer function from Control System Toolbox (default: 2nd-order filter)
%   t_span    - [t_start, t_end] evaluation interval (default: [0, 12.0])
%
% Outputs:
%   fourier_meta - Struct containing coefficients, w0, and Leibniz evaluation
%   sol_analytic - Struct containing function handle of truncated synthesis
%   sim_data     - Struct containing lsim simulation data [t, u_synth, y_out]

    arguments
        L (1,1) double {mustBePositive} = 2.0
        t_shift (1,1) double {mustBeNumeric} = 0.0
        Amp (1,1) double {mustBeNumeric} = 1.0
        N_terms (1,1) double {mustBeInteger, mustBePositive} = 15
        sys_plant = tf(9.0, [1.0, 1.5, 9.0]) % Default 2nd-order plant: wn=3, zeta=0.25
        t_span (1,2) double {mustBeNumeric} = [0, 12.0]
    end

    %% 1. Analytical Fourier Coefficient Generation
    % Fundamental circular frequency: w0 = pi / L
    w0 = pi / L;
    
    k_indices = 0:(N_terms - 1);
    harmonic_orders = 2 * k_indices + 1; % Odd harmonics: 1, 3, 5, ...
    
    % Unshifted square wave coefficients: bn_base = 4*Amp / (n*pi)
    bn_base = (4 * Amp) ./ (harmonic_orders * pi);
    
    % Incorporate time-shift: sin(n*w0*(t + t_shift)) = sin(n*w0*t)*cos(n*w0*t_shift) + cos(n*w0*t)*sin(n*w0*t_shift)
    % an = bn_base * sin(n*w0*t_shift)
    % bn = bn_base * cos(n*w0*t_shift)
    phase_shifts = harmonic_orders * w0 * t_shift;
    an_coeffs = bn_base .* sin(phase_shifts);
    bn_coeffs = bn_base .* cos(phase_shifts);

    fourier_meta.L = L;
    fourier_meta.Period = 2 * L;
    fourier_meta.w0 = w0;
    fourier_meta.harmonic_orders = harmonic_orders;
    fourier_meta.an = an_coeffs;
    fourier_meta.bn = bn_coeffs;

    % Leibniz Series Evaluation (Special case when L = pi, t_shift = 0)
    leibniz_sum = sum(((-1).^k_indices) ./ harmonic_orders);
    fourier_meta.leibniz_pi_estimate = 4 * leibniz_sum;

    %% 2. Function Handle Synthesis
    synth_handle = @(t) evaluate_general_fourier(t, w0, harmonic_orders, an_coeffs, bn_coeffs);
    sol_analytic.f_synth = synth_handle;

    %% 3. Numerical Simulation via Control System Toolbox (lsim)
    t_eval = linspace(t_span(1), t_span(2), 1200);
    u_eval = synth_handle(t_eval);

    if isa(sys_plant, 'tf') || isa(sys_plant, 'zpk') || isa(sys_plant, 'ss')
        [y_sim, ~] = lsim(sys_plant, u_eval, t_eval);
    else
        y_sim = u_eval;
    end

    sim_data.t = t_eval;
    sim_data.u = u_eval;
    sim_data.y = y_sim;

    %% 4. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Fourier Series II Solver: MIT 18.03 Lecture 21   \n');
    fprintf('====================================================\n');
    fprintf('Window Half-Period (L)  : %.4f s (Full Period = %.4f s)\n', L, 2*L);
    fprintf('Fundamental Freq (w0)   : %.4f rad/s\n', w0);
    fprintf('Time-Shift Offset       : %.4f s\n', t_shift);
    fprintf('Harmonic Truncation     : %d odd terms (up to harmonic %d)\n', ...
            N_terms, harmonic_orders(end));
    fprintf('Leibniz Pi Estimate     : %.6f (Error: %.2e)\n', ...
            fourier_meta.leibniz_pi_estimate, abs(fourier_meta.leibniz_pi_estimate - pi));
    fprintf('Plant Model G(s)        :\n');
    display(sys_plant);
    fprintf('====================================================\n\n');

    %% 5. Standalone Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 1050, 480]);

        % Left Subplot: Harmonic Amplitude Distribution
        subplot(1, 2, 1);
        c_mag = sqrt(an_coeffs.^2 + bn_coeffs.^2);
        stem(harmonic_orders, c_mag, 'filled', 'LineWidth', 1.5, 'MarkerSize', 5);
        grid on; box on;
        xlabel('Harmonic Order n', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Harmonic Magnitude \surd(a_n^2 + b_n^2)', 'FontSize', 10, 'FontWeight', 'bold');
        title('Harmonic Spectrum (1/n Envelope Decay)', 'FontSize', 11, 'FontWeight', 'bold');

        % Right Subplot: Time-Domain Dynamic Response via lsim
        subplot(1, 2, 2);
        hold on; grid on; box on;

        plot(t_eval, u_eval, 'r--', 'LineWidth', 1.4, ...
             'DisplayName', sprintf('Synthesized Input u(t) (N = %d)', N_terms));
        plot(t_eval, y_sim,  'b-',  'LineWidth', 2.0, ...
             'DisplayName', 'LTI Filtered Output y(t) = G(s)u(t)');

        xlabel('Time t (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('State Amplitude', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('Dynamic Tracking: Period = %.1fs, Shift = %.2fs', 2*L, t_shift), ...
              'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 8.5);
        xlim([t_span(1), t_span(2)]);
        hold off;
    end
end

%% Local Helper for Fourier Evaluation
function val = evaluate_general_fourier(t_in, w0, orders, an, bn)
    val = zeros(size(t_in));
    for idx = 1:length(orders)
        n = orders(idx);
        val = val + an(idx) * cos(n * w0 * t_in) + bn(idx) * sin(n * w0 * t_in);
    end
end