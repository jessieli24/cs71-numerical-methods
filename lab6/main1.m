% 
% main1.m - solves the Lotka-Volterra equations with RK4, RK2, and 
%       Euler's method
% 
% Jessie Li, CS 71 Fall 2023
%

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

% -------------------- constants ---------------------- %
% ----------------------------------------------------- %
r = 1.1;    % intrinsic rate of population increase of prey
a = 0.05;   % predation rate coefficient
b = 0.01;   % reproduction rate of predators born per prey eaten
m = 0.4;    % mortality rate of the predators

tMin = 0;
tMax = 50;

dH_dt = @(H, F, t) r * H - a * H * F;
dF_dt = @(H, F, t) b * H * F - m * F;

f = @(y, t) [dH_dt(y(1), y(2), t); 
             dF_dt(y(1), y(2), t)];

H0 = 250;   % initial prey population
F0 = 10;    % initial predator population

x0 = [H0; F0];

minSteps = 4;
maxSteps = 12;

% number of step sizes to try
nSteps = maxSteps - minSteps + 1;

% ----------------------------------------------------- %
% H(t) and F(t) for t = 0 to 50
% ----------------------------------------------------- %
% step sizes, h = 2^-n
h = zeros(1, nSteps);

% store estimate at t = 50 for each step size
w50_rk4 = zeros(2, nSteps);
w50_rk2 = zeros(2, nSteps);
w50_euler = zeros(2, nSteps);

for n = minSteps : maxSteps
    i = n - minSteps + 1; % index

    h(i) = 2^-n;

    % calculate H(50) and F(50)
    w_rk4 = rk4(f, tMin, tMax, h(i), x0);
    w50_rk4(:, i) = w_rk4(:, end);

    w_rk2 = rk2(f, tMin, tMax, h(i), x0);
    w50_rk2(:, i) = w_rk2(:, end);

    w_euler = euler(f, tMin, tMax, h(i), x0);
    w50_euler(:, i) = w_euler(:, end);

    % graph the solution only for the most accurate and least accurate step sizes
    if n == minSteps || n == maxSteps
        plotLVSolution(w_rk4, w_rk2, w_euler, tMin, tMax, n)
    end
end

% ----------------------------------------------------- %
% plot parametric solution < H(t), F(t) >
% ----------------------------------------------------- %
figure 
setColors()

hold on
plot(w_rk4(1, :), w_rk4(2, :))
plot(w_rk2(1, :), w_rk2(2, :))
plot(w_euler(1, :), w_euler(2, :))
hold off

xlabel('H')
ylabel('F')

legend('RK4', 'RK2', 'Euler')
% ----------------------------------------------------- %
% plot error at t = 50 v. step size
% ----------------------------------------------------- %
% assume the exact solution is RK4 with the smallest step size (h = 2^-12)
y50_H = w50_rk4(1, end);
y50_F = w50_rk4(2, end);

errors_rk4_H = abs(y50_H - w50_rk4(1, :));
errors_rk4_F = abs(y50_F - w50_rk4(2, :));

errors_rk2_H = abs(y50_H - w50_rk2(1, :));
errors_rk2_F = abs(y50_F - w50_rk2(2, :));

errors_euler_H = abs(y50_H - w50_euler(1, :));
errors_euler_F = abs(y50_F - w50_euler(2, :));

figure
setColors()

hold on
plot(minSteps : maxSteps, errors_rk4_H, '-o')
plot(minSteps : maxSteps, errors_rk4_F, '-o')

plot(minSteps : maxSteps, errors_rk2_H, '-o')
plot(minSteps : maxSteps, errors_rk2_F, '-o')

plot(minSteps : maxSteps, errors_euler_H, '-o')
plot(minSteps : maxSteps, errors_euler_F, '-o')

hold off

set(gca, 'YScale', 'log')
xlabel('n (step size = $2^{-n}$ )')
ylabel('Absolute error')
title({'Absolute Error of H(t) and F(t)', 'v. Step Size at t = 50'})
legend('RK4, H', 'RK4, F', 'RK2, H', 'RK2, F',  'Euler, H', 'Euler, F')

% ---------------- helper functions ------------------- %
% ----------------------------------------------------- %

function plotLVSolution(w_rk4, w_rk2, w_euler, a, b, n)
    h = 2^-n;

    figure 
    setColors()

    hold on
    plot(a:h:b, w_rk4(1, :), 'LineWidth', 3)
    plot(a:h:b, w_rk4(2, :), 'LineWidth', 3)

    plot(a:h:b, w_rk2(1, :), '--', 'LineWidth', 2)
    plot(a:h:b, w_rk2(2, :), '--', 'LineWidth', 2)

    plot(a:h:b, w_euler(1, :), '-.', 'LineWidth', 2)
    plot(a:h:b, w_euler(2, :), '-.', 'LineWidth', 2)
    hold off

    xlabel('t')
    legend('RK4, H', 'RK4, F', 'RK2, H', 'RK2, F',  'Euler, H', 'Euler, F')

     title(sprintf('H(t) and F(t) (step size = $2^{-%d} )$', n))
end

function setColors()
    color_order = [0.37 0.60 0.94
                   0.05 0.26 0.57
                   0.98 0.58 0.89
                   0.99 0.82 0.54
                   0.81 0.59 0.95
                   0.53 0.98 0.84];
    
    colororder(color_order)
end


