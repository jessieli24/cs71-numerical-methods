% 
% main2.m - parametric polynomial interpolation with linearly-spaced points
% 
% Jessie Li, CS 71 Fall 2023
%

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

data = readmatrix('lab3_data.txt');

sdata = data(:, 1);
xdata = data(:, 2);
ydata = data(:, 3);

for ds = 1 : 6
    s = sdata(1:ds:end);
    x = xdata(1:ds:end);
    y = ydata(1:ds:end);
    
    a = min(s);     % left bound of interval
    b = max(s);     % right bound of interval
    h = 0.001;      % step size for plotting

    ss = a:h:b;
    
    xx = arrayfun(@(sp) neville(s, x, length(s), sp), ss);
    yy = arrayfun(@(sp) neville(s, y, length(s), sp), ss);

    fprintf('------ ds = %d, (n = %d) ------\n', ds, length(s) - 1);
    fprintf('min x: %0.5f, max x: %0.5f\n', min(xx), max(xx));
    fprintf('min y: %0.5f, max y: %0.5f\n', min(yy), max(yy));
    fprintf('-------------------------------\n');
    
    figure
    t = tiledlayout(2, 1);

    % plot x(s)
    nexttile
    plot(s, x, 'bo', ss, xx, 'b-', sdata, xdata, 'kx')
    ylabel('x(s)')

    % plot y(s)
    nexttile
    plot(s, y, 'mo', ss, yy, 'm-', sdata, ydata, 'kx')
    ylabel('y(s)')

    % label
    title(t, sprintf('$x(s)$ and $y(s)$ for $\\Delta s = %d$ (Lagrange)', ds));
    xlabel(t, 's')
    
    % saveas(gcf, sprintf('results/q2-ds%d-sxsy.png', ds));
    
    % plot interpolated shape
    figure
    plot(x, y, 'o', xx, yy, '-', xdata, ydata, 'kx')
    
    xlabel('x')
    ylabel('y')
    title(sprintf('Lagrange Interpolation ($\\Delta s = %d$)', ds));
    
    % saveas(gcf, sprintf('results/q2-ds%d-xy.png', ds));
end
