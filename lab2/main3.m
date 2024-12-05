%
% main3.m - estimates the root of f for Q3B
% 
% Jessie Li, CS 71 Fall 2023
%

function main3(method)

% set default font to Times New Roman for all graphs
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

% -------------------- constants ---------------------- %
% ----------------------------------------------------- %
filename = 'results/q3';
graphTitle = 'Nonlinear Newton';

MAX_ITER = 100;

% TOL = 2.2204e-16
% MATLAB's default for root-finding
% https://www.mathworks.com/help/matlab/math/setting-options.html#bt00l89-1
TOL = optimset('fzero').TolX;
% ----------------------------------------------------- %
% ----------------------------------------------------- %

x0 = [deg2rad(59); deg2rad(348)];

[n, x, err] = newtonNonlinear(@f, @J, x0, MAX_ITER, TOL);

% plot errors v. iterations
figure
semilogy(1:n, err(1:n));

ylabel('Relative error');
xlabel('Number of iterations');
title(graphTitle);

% iteration number must be an integer
xlabels = get(gca, 'xTick');
xticks(unique(round(xlabels)));

saveas(gcf, sprintf('%s.png', filename));

% compute rate of convergence (expect alpha = 2)
a = calculateConvergence(n, err);

figure
plot(3:n, a);

ylabel('Rate of convergence');
xlabel('Number of iterations');
title(graphTitle);

% iteration number must be an integer
xlabels = get(gca, 'xTick');
xticks(unique(round(xlabels)));

saveas(gcf, sprintf('%s-conv.png', filename));
end

function y = f(x)
r1 = 45;
r2 = 32;
r3 = 33;
r4 = 21;

theta1 = deg2rad(80);
theta2 = x(1);
theta3 = x(2);
theta4 = theta1 + pi;

y = [r2*cos(theta2) + r3*cos(theta3) + r4*cos(theta4) - r1;
    r2*sin(theta2) + r3*sin(theta3) + r4*sin(theta4)];
end

function y = J(x)
r2 = 32;
r3 = 33;

theta2 = x(1);
theta3 = x(2);

y = [-r2*sin(theta2), -r3*sin(theta3);
    r2*cos(theta2),  r3*cos(theta3)];
end
