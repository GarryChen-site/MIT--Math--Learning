%% Test Driver Script: MIT 18.03 Lecture 1 Question 1
clear; clc; close all;

% Problem Parameters: dy/dx = 2y + 1, y(0) = 0
a = 2;
b = 1;
x0 = 0;
y0 = 0;
x_span = [0, 1.5];

% Solve using the modular solver
[y_exact, num_sol] = solve_first_order_affine_ode(a, b, x0, y0, x_span);

% Plot comparison
figure('Color', [1 1 1]);
x_fine = linspace(x_span(1), x_span(2), 200);
plot(x_fine, y_exact(x_fine), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytic Solution');
hold on; grid on;
plot(num_sol.x, num_sol.y, 'ro', 'MarkerSize', 5, 'DisplayName', 'Numeric ode45');

xlabel('x', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('y(x)', 'FontSize', 11, 'FontWeight', 'bold');
title("Solution Comparison: dy/dx = 2y + 1, y(0) = 0", 'FontSize', 12);
legend('Location', 'northwest');