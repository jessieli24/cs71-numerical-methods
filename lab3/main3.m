% 
% main3.m - parametric cubic spline interpolation
% 
% Jessie Li, CS 71 Fall 2023
%

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

data = readmatrix('lab3_data.txt');

s = data(:, 1);
x = data(:, 2);
y = data(:, 3);

a = min(s);     % left bound of interval
b = max(s);     % right bound of interval
n = length(s);  % number of points
m = 1000;       % how many points to plot per segment

[~, xx] = evaluateSpline(n, s, x, m, 'results/q3-ds1-x.csv');
[ss, yy] = evaluateSpline(n, s, y, m, 'results/q3-ds1-y.csv');

fprintf('------ ds = 1 ------\n');
fprintf('min x: %0.5f, max x: %0.5f\n', min(min(xx)), max(max(xx)));
fprintf('min y: %0.5f, max y: %0.5f\n', min(min(yy)), max(max(yy)));

figure
t = tiledlayout(2, 1);

% plot x(s)
nexttile
plot(ss, xx, 'b-', s, x, 'bo')
ylabel('x(s)')

% plot y(s)
nexttile
plot(ss, yy, 'm-', s, y, 'mo')
ylabel('y(s)')

title(t, 'Cubic Spline Interpolation of $x(s)$ and $y(s)$', 'interpreter', 'latex');
xlabel(t, 's', 'interpreter', 'latex');

% saveas(gcf, 'results/q3-ds1-sxsy.png');

% plot interpolated shape
figure
plot(xx, yy, '-', x, y, 'o')
   
xlabel('x')
ylabel('y')
title('Cubic Spline Interpolation ($\Delta s = 1$)');

% saveas(gcf, 'results/q3-ds1-xy.png');

% Returns a cubic spline interpolation of the x, y points. 
% 
%     Parameters:
%         n: number of data points
%         x: x values as a list of size n
%         y: y values as a list of size n
%         m: number of points to plot per segment
%         filename: name of csv for saving coefficients 

function [xx, yy] = evaluateSpline(n, x, y, m, filename)   
    % get coefficients for each segment 
    [a, b, c, d] = cubicSpline(n, x, y);
    % writematrix([a b c d], filename);
    
    % store interpolated points in xx and yy
    xx = zeros(m, n-1);
    yy = zeros(m, n-1);
    
    % get interpolated points for each segment
    for i = 1 : n-1
        
        % cubic interpolation for this segment
        p = @(xp) a(i) ...
                + b(i).*(xp - x(i)) ...
                + c(i).*(xp - x(i)).^2 ...
                + d(i).*(xp - x(i)).^3;
        
        xx(:, i) = linspace(x(i), x(i+1), m);
        yy(:, i) = p(xx(:, i)); 
    end  
end
