% 
% main2.m - compares the stability of the solution to y' = -7y 
%           for several methods by varying step size
%
% Jessie Li, CS 71 Fall 2023
%

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

% -------------------- constants ---------------------- %
% ----------------------------------------------------- %
f = @(y, t) -7 * y;

% analytical solution
y0 = 50;
y = @(t) y0 .* exp(-7 .* t);

tMin = 0;
tMax = 50;

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% calculated maximum allowable step-size for each method:
%   h = 1/7 (AB2)
%   h = 2/7 (RK2)
%   h = 12/35 (AB2/AM2)
%   h = 6/7 (AM2)

for h = [1/8, 1/7, 1/6, 1/3, 1/2, 1]
    w_rk2 = rk2(f, tMin, tMax, h, y0);
    w_ab2 = ab2(f, tMin, tMax, h, y0, y);
    w_am2 = am2Particular(tMin, tMax, h, y0, y);
    w_adams2 = adams2(f, tMin, tMax, h, y0, y);

    figure
    setColors()
    hold on
    
    if h < 1/3
        plot(tMin:1e-3:tMax, y(tMin:1e-3:tMax), 'LineWidth', 2)
        plot(tMin:h:tMax, w_am2, '--')
        plot(tMin:h:tMax, w_adams2, '--')
        plot(tMin:h:tMax, w_rk2, '--')
        plot(tMin:h:tMax, w_ab2, '--')
    
        legend('Exact', 'AM2', 'AB2/AM2', 'RK2', 'AB2')
        xlim([0 3])

    elseif h < 1/2
        plot(tMin:1e-3:tMax, y(tMin:1e-3:tMax), 'LineWidth', 2)
        plot(tMin:h:tMax, w_am2, '--')
        plot(tMin:h:tMax, w_adams2, '--')
        plot(tMin:h:tMax, w_rk2, '--')
    
        legend('Exact', 'AM2', 'AB2/AM2', 'RK2')
        xlim([0 3])
    
    elseif h < 1
        plot(tMin:1e-3:tMax, y(tMin:1e-3:tMax), 'LineWidth', 2)
        plot(tMin:h:tMax, w_am2, '--')
        plot(tMin:h:tMax, w_adams2, '--')
    
        legend('Exact', 'AM2', 'AB2/AM2')
        xlim([0 2])
    
    elseif h < 2
        plot(tMin:1e-3:tMax, y(tMin:1e-3:tMax), 'LineWidth', 2)
        plot(tMin:h:tMax, w_am2, '--')
    
        legend('Exact', 'AM2')
        xlim([0 10])
    end
    
    hold off

    xlabel('t')
    ylabel('y')
    title(sprintf('h = 1/%d', 1/h))
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
