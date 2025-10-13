%% YourLastName-Lab1-main.m
% Lab 1 - Differential Equations with MATLAB
% Name: Nouh Shaikh
% Date: October 12 2025
% Course: MTH 286
%
% ------------------------------------------------------------
% Problem 1: Solve dy/dx = sin(x) - 2y, y(0) = 2
% ------------------------------------------------------------

clc; clear; close all

%% (a) Exact Solution Derivation (symbolic form)
% dy/dx + 2y = sin(x)
% Integrating factor = e^(2x)
% Exact solution: y(x) = (2*sin(x) - cos(x))/5 + (11/5)*exp(-2x)

y_exact = @(x) (2*sin(x) - cos(x))/5 + (11/5)*exp(-2*x);

%% (b) Evaluate y(x) at x = 1, 5, 15
x_vals = [1 5 15];
y_vals = y_exact(x_vals);

fprintf('Problem 1(b): Exact solution values\n');
for i = 1:length(x_vals)
    fprintf('y(%2.0f) = %.6f\n', x_vals(i), y_vals(i));
end

%% (c) Plot exact and numerical solutions on [0, 25]
f = @(x, y) sin(x) - 2*y;  
[x_num, y_num] = ode45(f, [0 25], 2); 

x_exact = linspace(0, 25, 200);
y_exact_vals = y_exact(x_exact);

figure(1)
plot(x_exact, y_exact_vals, 'b-', 'LineWidth', 1.5); hold on
plot(x_num, y_num, 'ro', 'MarkerSize', 3);
xlabel('x'); ylabel('y(x)');
title('Problem 1(c): Exact vs Numerical Solution of dy/dx = sin(x) - 2y');
legend('Exact Solution', 'Numerical (ode45)', 'Location', 'northeast');
grid on

%% ------------------------------------------------------------
% Problem 2: Logistic Equation
% dP/dt = P(1 - P),   0 ≤ t ≤ 4,   0 ≤ P ≤ 5
% ------------------------------------------------------------

[t, P] = meshgrid(0:0.2:4, 0:0.2:5); 
dPdt = P .* (1 - P);
dt = ones(size(t));

figure(2)
quiver(t, P, dt, dPdt, 'r');
xlabel('t'); ylabel('P');
title('Problem 2: Slope Field for Logistic Equation dP/dt = P(1 - P)');
hold on

% Plot three solutions with different initial conditions
f2 = @(t, P) P .* (1 - P);
P0 = [0.2, 0.7, 5];
colors = ['b', 'g', 'k'];

for i = 1:length(P0)
    [t_num, P_num] = ode45(f2, [0 4], P0(i));
    plot(t_num, P_num, 'Color', colors(i), 'LineWidth', 1.5);
end

legend('Slope Field', 'P(0)=0.2', 'P(0)=0.7', 'P(0)=5', 'Location', 'east');
grid on
hold off

fprintf('\nProblem 2 completed: Slope field and trajectories plotted.\n');
