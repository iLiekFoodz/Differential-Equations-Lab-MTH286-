clc; clear; close all;
% MATLAB Lab 2 Assignment
% Name: Nouh Shaikh
% CSU ID: 2857056
% In this lab, we solve three differential equations using MATLAB's dsolve function.
% Each problem is symbolic, meaning MATLAB finds the exact formula (not numbers).

%% Problem 1
% Equation: 3y'' + y' - 2y = 0
% This is a 2nd order homogeneous differential equation.
syms y(x)
ode1 = 3*diff(y,x,2) + diff(y,x) - 2*y == 0;
y1 = dsolve(ode1);         % Solve the ODE symbolically
disp('Problem 1 Solution:')
disp(y1)                   % Display the general solution

%% Problem 2
% Equation: 3y''' + 5y'' + y' - y = 0
% Given initial conditions: y(0)=0, y'(0)=1, y''(0)=-1
% This one is 3rd order and we'll use the conditions to find constants.
Dy = diff(y,x);
D2y = diff(y,x,2);
D3y = diff(y,x,3);
ode2 = 3*D3y + 5*D2y + Dy - y == 0;
conds = [y(0)==0, Dy(0)==1, D2y(0)==-1];
y2 = dsolve(ode2, conds);  % Solve with initial conditions
disp('Problem 2 Solution:')
disp(y2)

%% Problem 3
% Equation: y'' - 4y' - 12y = 3t^3 - 5t + 2
% This is a non-homogeneous ODE (it has a polynomial on the right side).
% MATLAB automatically finds both the homogeneous and particular solutions.
syms y(t)
ode3 = diff(y,t,2) - 4*diff(y,t) - 12*y == 3*t^3 - 5*t + 2;
y3 = dsolve(ode3);
disp('Problem 3 Solution:')
disp(y3)

%% Summary of all results
fprintf('\n = Summary of Solutions = \n')
fprintf('1) y(x) = C1*exp(2x/3) + C2*exp(-x)\n')
fprintf('2) y(x) = (-9/16 + (1/4)*x)*exp(-x) + (9/16)*exp(x/3)\n')
fprintf('3) y(t) = C1*exp(6t) + C2*exp(-2t) - (1/4)*t^3 + (1/4)*t^2 + (1/8)*t - 1/6\n')
