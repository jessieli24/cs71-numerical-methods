% 
% main1.m - interpolates polynomials with linearly-spaced and Chebyshev
% nodes
% 
% Jessie Li, CS 71 Fall 2023
%

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

% visualizing the function 

a = -2;         % left bound of interval
b = 2;          % right bound of interval
h = 0.001;      % step size for plotting

xx = a:h:b;
yy = f(xx);

plot(xx, yy)

xlabel('x');
ylabel('f(x)');

% interpolation 

for n = 5:5:30 % n = order of polynomial
    % linear -- input points
    xl = linspace(a, b, n+1);
    yl = f(xl);

    % linear -- interpolated points
    yyl = arrayfun(@(x) neville(xl, yl, n+1, x), xx);

    % Chebyshev -- input points
    xc = chebyshevPoints(a, b, n+1);
    yc = f(xc);

    % Chebyshev -- plotting points
    yyc = arrayfun(@(x) neville(xc, yc, n+1, x), xx);

    % plot interpolations
    figure

    hold on
    plot(xx, yy);
    plot(xl, yl, 'go', xx, yyl, 'g--');
    plot(xc, yc, 'mo', xx, yyc, 'm--');
    hold off
    
    legend('f', '', 'linear', '', 'chebyshev', 'Location', 'southeast');

    t = sprintf('Polynomial Interpolations of $f$ (n = %d)', n);
    title(t);

    xlabel('x');
    ylabel('y');
    
    % saveas(gcf, sprintf('results/q1-interp-%d.png', n));
    
    % plot errors
    figure

    hold on
    plot(xx, abs(yy - yyl), 'g');
    plot(xx, abs(yy - yyc), 'm');
    hold off

    legend('linear', 'chebyshev', 'Location', 'southeast');

    t = sprintf('Interpolation Errors (n = %d)', n);
    title(t);

    xlabel('x');
    ylabel('Absolute error');
    
    set(gca, 'YScale', 'log')
    
    % saveas(gcf, sprintf('q1-err-%d.png', n));
    
    fprintf('------------ n = %d ------------\n', n);
    fprintf('max linear error: %0.5f\n', max(abs(yy-yyl)));
    fprintf('max chebyshev error: %0.5f\n', max(abs(yy-yyc)));
    fprintf('--------------------------------\n');
end

function y = f(x)
    y = sin(6.*x) .* cos(5.*x) - x.^2 .* exp(-x./5);
end
