function [result_val, roots_array, polar_info] = solve_complex_power_roots(z_input, n_degree, calc_mode)
% SOLVE_COMPLEX_POWER_ROOTS Evaluates powers and n-th roots of complex numbers.
%
% % Required Toolboxes: 
% %   - Base MATLAB (Zero external dependency).
%
% % MBD Concept Mapping: Demonstrates frequency response phasor transformations,
% %                      multiplication of magnitudes, addition of angles, and
% %                      distribution of poles along radial geometries in the s/z planes.
%
% Inputs (All optional with defaults matching MIT 18.03 Lec 5 Question 3):
%   z_input   - Input complex number (default: 1 + 1i)
%   n_degree  - Integer power or root order (default: 4)
%   calc_mode - String: 'power' computes z^n; 'roots' computes all n roots of z (default: 'power')
%
% Outputs:
%   result_val  - Complex result for 'power' mode (NaN for 'roots' mode)
%   roots_array - Array of all n distinct complex roots for 'roots' mode
%   polar_info  - Struct containing magnitude and phase information

    arguments
        z_input (1,1) double = 1 + 1i
        n_degree (1,1) double {mustBeInteger, mustBePositive} = 4
        calc_mode (1,1) string {mustBeMember(calc_mode, ["power", "roots"])} = "power"
    end

    %% 1. Polar Representation Decomposition
    r_modulus = abs(z_input);
    theta_rad = angle(z_input);
    theta_deg = theta_rad * (180 / pi);

    polar_info.input_modulus = r_modulus;
    polar_info.input_phase_rad = theta_rad;
    polar_info.input_phase_deg = theta_deg;

    %% 2. Power and Roots Computation via De Moivre's Theorem
    if calc_mode == "power"
        % Powers: |z^n| = |z|^n, Arg(z^n) = n * Arg(z)
        pow_modulus = r_modulus^n_degree;
        pow_phase_rad = n_degree * theta_rad;
        
        result_val = pow_modulus * exp(1i * pow_phase_rad);
        roots_array = [];
        
        polar_info.result_modulus = pow_modulus;
        polar_info.result_phase_rad = angle(result_val); % Wrapped to [-pi, pi]
        polar_info.result_phase_deg = polar_info.result_phase_rad * (180 / pi);
        
    else % 'roots' mode: w_k = r^(1/n) * exp(j * (theta + 2*pi*k)/n)
        result_val = NaN;
        root_modulus = r_modulus^(1 / n_degree);
        roots_array = zeros(1, n_degree);
        
        for k = 0:(n_degree - 1)
            phi_k = (theta_rad + 2 * pi * k) / n_degree;
            roots_array(k + 1) = root_modulus * exp(1i * phi_k);
        end
        
        polar_info.result_modulus = root_modulus;
        polar_info.result_phase_rad = angle(roots_array);
        polar_info.result_phase_deg = polar_info.result_phase_rad * (180 / pi);
    end

    %% 3. Terminal Diagnostics Output
    fprintf('====================================================\n');
    fprintf('   Complex Number Operator: De Moivre Solver        \n');
    fprintf('====================================================\n');
    fprintf('Input Complex z        : %.4f %+0.4fi\n', real(z_input), imag(z_input));
    fprintf('Modulus |z|            : %.4f\n', r_modulus);
    fprintf('Argument Arg(z)        : %.4f rad (%.2f deg)\n', theta_rad, theta_deg);
    fprintf('Degree / Order n       : %d\n', n_degree);
    fprintf('Operating Mode         : %s\n', calc_mode);
    fprintf('----------------------------------------------------\n');
    
    if calc_mode == "power"
        fprintf('Result z^%d            : %.4f %+0.4fi\n', n_degree, real(result_val), imag(result_val));
        fprintf('Result Modulus         : %.4f\n', polar_info.result_modulus);
        fprintf('Result Phase (Wrapped) : %.4f rad (%.2f deg)\n', ...
                polar_info.result_phase_rad, polar_info.result_phase_deg);
    else
        fprintf('Evaluated %d Distinct Roots:\n', n_degree);
        for idx = 1:n_degree
            fprintf('  Root %d: %8.4f %+0.4fi  | Phase: %6.1f deg\n', ...
                    idx, real(roots_array(idx)), imag(roots_array(idx)), ...
                    polar_info.result_phase_deg(idx));
        end
    end
    fprintf('====================================================\n\n');

    %% 4. Demonstration Visualizer (Active on direct F5 run)
    if nargout == 0
        figure('Color', [1 1 1], 'Position', [150, 150, 800, 520]);
        hold on; grid on; box on;
        
        xline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        if calc_mode == "power"
            % Plot input vector
            quiver(0, 0, real(z_input), imag(z_input), 0, 'b-', 'LineWidth', 1.8, ...
                   'DisplayName', sprintf('Input z = %.1f+%.1fi', real(z_input), imag(z_input)));
            % Plot power vector
            quiver(0, 0, real(result_val), imag(result_val), 0, 'r-', 'LineWidth', 2.0, ...
                   'DisplayName', sprintf('Result z^{%d} = %.1f+%.1fi', n_degree, real(result_val), imag(result_val)));
            
            plot(real(result_val), imag(result_val), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
            title(sprintf('De Moivre Power: (%.1f + %.1fi)^{%d} = %.1f %+0.1fi', ...
                  real(z_input), imag(z_input), n_degree, real(result_val), imag(result_val)), ...
                  'FontSize', 11, 'FontWeight', 'bold');
        else
            % Plot circular boundary for root distribution
            rad = polar_info.result_modulus;
            circ_theta = linspace(0, 2*pi, 200);
            plot(rad * cos(circ_theta), rad * sin(circ_theta), 'k:', 'LineWidth', 1.2, ...
                 'DisplayName', sprintf('Root Circle (R = %.2f)', rad));
            
            % Plot root vectors
            for idx = 1:n_degree
                quiver(0, 0, real(roots_array(idx)), imag(roots_array(idx)), 0, ...
                       'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
            plot(real(roots_array), imag(roots_array), 's', 'MarkerEdgeColor', 'r', ...
                 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'DisplayName', 'Computed Roots');
             
            title(sprintf('%d Roots of Complex Number z = %.1f + %.1fi', ...
                  n_degree, real(z_input), imag(z_input)), 'FontSize', 11, 'FontWeight', 'bold');
        end

        axis equal;
        xlabel('Real Axis Re(z)', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Imaginary Axis Im(z)', 'FontSize', 10, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        hold off;
    end
end