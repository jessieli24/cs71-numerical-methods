%
% main2b.m - estimates the root of f for Q2B using the specified method 
% 
% Input:
%   method (string):
%       * 'newton'
%       * 'secant'
%       * 'newtonModified'
%       * 'newtonCubic'
% 
% Jessie Li, CS 71 Fall 2023
%

function main2b(method)

% set default font to Times New Roman for all graphs
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

% -------------------- constants ---------------------- %
% ----------------------------------------------------- %
filename = sprintf('results/q2b-%s', method);

% interval
a = 1;
b = 2;

MAX_ITER = 100;

% TOL = 2.2204e-16
% MATLAB's default for root-finding
% https://www.mathworks.com/help/matlab/math/setting-options.html#bt00l89-1
TOL = optimset('fzero').TolX;
% ----------------------------------------------------- %
% ----------------------------------------------------- %

x0 = (a + b)/2;

switch method
    case 'newton'
        graphTitle = 'Newton';
        [n, x, err] = newton(@f, @df, x0, MAX_ITER, TOL);
    case 'secant'
        graphTitle = 'Secant';
        [n, x, err] = secant(@f, a, b, MAX_ITER, TOL);
    case 'newtonModified'
        graphTitle = 'Modified Newton';
        [n, x, err] = newtonModified(@f, @df, @ddf, x0, MAX_ITER, TOL);
    case 'newtonCubic'
        graphTitle = 'Cubic Newton';
        [n, x, err] = newtonCubic(@f, @df, @ddf, x0, MAX_ITER, TOL);
    otherwise
        error('Unknown method. Must be `newton`, `secant`, `newtonModified`, or `newtonCubic`.');
end

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

% compute rate of convergence
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
y = ((x + cos(x)) * exp(-x^2) + x*cos(x))^2;
end

function y = df(x)
y = -2 * (exp(-x^2) * (x + cos(x)) + x * cos(x)) * ...
    (exp(-x^2) * (sin(x) - 1) - cos(x) + x * sin(x) + ...
    2 * x * exp(-x^2) * (x + cos(x)));
end

function y = ddf(x)
y = 2 * (exp(-x^2) * (sin(x) - 1) - cos(x) + x * sin(x) + ...
    2 * x * exp(-x^2) * (x + cos(x)))^2 - 2 * (exp(-x^2) * ...
    (x + cos(x)) + x * cos(x)) * (2 * sin(x) + 2 * exp(-x^2) * ...
    (x + cos(x)) + x * cos(x) + exp(-x^2) * cos(x) - 4 * x^2 * ...
    exp(-x^2) * (x + cos(x)) - 4 * x * exp(-x^2) * (sin(x) - 1));
end
