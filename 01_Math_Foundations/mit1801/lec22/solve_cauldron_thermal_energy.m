function sol = solve_cauldron_thermal_energy(a_height_m, T_bottom_degC, T_top_degC, candy_bar_kcal)
% SOLVE_CAULDRON_THERMAL_ENERGY Calculates the enclosed liquid volume of a paraboloid
% cauldron (y = x^2 up to height a) and the thermal energy required to heat the liquid
% given a linear temperature gradient T(y) from bottom to top.
% Formulas:
%   Volume       : V = pi * a^2 / 2  (m^3)
%   Thermal Work : H = \int_0^a T(y) * (pi * y) dy  (degC * m^3)
%   Calories     : Cal = H * (100 cm / 1 m)^3 * (1 cal / (cm^3 * degC))
%
% Syntax:
%   sol = solve_cauldron_thermal_energy(a_height_m, T_bottom_degC, T_top_degC, candy_bar_kcal)
%
% Inputs:
%   a_height_m      - Maximum liquid height a in meters (default = 1.0 m)
%   T_bottom_degC   - Temperature at cauldron base y = 0 (default = 100 deg C)
%   T_top_degC      - Temperature at cauldron surface y = a (default = 70 deg C)
%   candy_bar_kcal  - Energy rating of reference candy bar (default = 250 kcal)
%
% Outputs:
%   sol             - Struct containing volume, thermal integral, calories, and candy bar equivalents

arguments
    a_height_m     (1,1) double {mustBePositive} = 1.0
    T_bottom_degC  (1,1) double = 100.0
    T_top_degC     (1,1) double = 70.0
    candy_bar_kcal (1,1) double {mustBePositive} = 250.0
end

% 1. Enclosed Volume via Disk/Shell Method: V = pi * a^2 / 2
volume_m3 = (pi * (a_height_m^2)) / 2.0;
volume_liters = volume_m3 * 1000.0;

% 2. Linear Temperature Gradient: T(y) = T_bottom - ( (T_bottom - T_top) / a ) * y
% Let slope = (T_bottom - T_top) / a
delta_T = T_bottom_degC - T_top_degC;
slope_T = delta_T / a_height_m;

% Thermal Integral: H = pi * \int_0^a [ T_bottom*y - slope_T * y^2 ] dy
% Antiderivative: pi * [ T_bottom * (y^2/2) - slope_T * (y^3/3) ] from 0 to a
H_val_m3_degC = pi * (T_bottom_degC * (a_height_m^2 / 2.0) - slope_T * (a_height_m^3 / 3.0));

% 3. Unit Conversions:
% 1 m^3 = 10^6 cm^3. Water specific heat capacity = 1 cal / (cm^3 * deg C)
total_calories = H_val_m3_degC * 1.0e6;
total_kcal = total_calories / 1000.0;
candy_bars_equivalent = total_kcal / candy_bar_kcal;

% Assemble solution structure
sol.liquid_height_m = a_height_m;
sol.volume_cubic_meters = volume_m3;
sol.volume_liters = volume_liters;
sol.bottom_temperature = T_bottom_degC;
sol.top_temperature = T_top_degC;
sol.thermal_integral_deg_m3 = H_val_m3_degC;
sol.total_calories = total_calories;
sol.total_kcal = total_kcal;
sol.equivalent_candy_bars = candy_bars_equivalent;

% Console report
fprintf('=== Cauldron Enclosed Volume & Thermal Energy Analysis ===\n');
fprintf('Cauldron Liquid Level (a)   : %.2f m\n', a_height_m);
fprintf('Total Liquid Volume (V)     : %.4f m^3 (%.1f Liters)\n', volume_m3, volume_liters);
fprintf('Temperature Gradient Field  : Bottom = %.1f deg C | Top = %.1f deg C\n', ...
    T_bottom_degC, T_top_degC);
fprintf('Thermal Integral (H)        : %.4f degC*m^3\n', H_val_m3_degC);
fprintf('Total Energy in Kilocalories: %.2e kcal\n', total_kcal);
fprintf('Equivalent Food Energy Units: ~%.1f Candy Bars (at %.0f kcal/bar)\n', ...
    candy_bars_equivalent, candy_bar_kcal);
end