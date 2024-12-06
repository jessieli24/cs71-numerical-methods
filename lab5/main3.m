% 
% main3.m - IVP with Euler's, Midpoint Rule (RK2), Modified Euler's, 
%           2-Step A-B/A-M Predictor-Corrector, and RK4
% 
% Jessie Li, CS 71 Fall 2023
%

% -------------------- constants ---------------------- %
% ----------------------------------------------------- %
% interval a <= t <= b
a = 0;
b = 3;

% initial value problem
dydt = @(y, t) -2 * y * t / (1 + t^2);
y0 = 1;

% analytical solution
y = @(t) 1 ./ (1 + t.^2);

n_min = 3;
n_max = 12;

% ---------------------- euler ------------------------ %
% ----------------------------------------------------- %
y2_error_euler = approximate(dydt, y0, y, a, b, n_min, n_max, 'euler');
% saveas(gcf, 'results/q3-euler.png');

% ------------------ midpoint/rk2 --------------------- %
% ----------------------------------------------------- %
y2_error_midpoint = approximate(dydt, y0, y, a, b, n_min, n_max, 'rk2');
% saveas(gcf, 'results/q3-midpoint.png')

% ---------------- modified euler --------------------- %
% ----------------------------------------------------- %
y2_error_modified_euler = approximate(dydt, y0, y, a, b, n_min, n_max, 'eulerModified');
% saveas(gcf, 'results/q3-modifiedEuler.png')

% ----------------- 2-step ab/am ---------------------- %
% ----------------------------------------------------- %
y2_error_adams = approximate(dydt, y0, y, a, b, n_min, n_max, 'adams2');
% saveas(gcf, 'results/q3-adams.png')

% ------------------------ rk4 ------------------------ %
% ----------------------------------------------------- %
y2_error_runge_kutta = approximate(dydt, y0, y, a, b, n_min, n_max, 'rk4');
% saveas(gcf, 'results/q3-rk4.png')

% --------------- all y(2) v. 1/∆t -------------------- %
% ----------------------------------------------------- %
inverse_delta_t = 2 .^ (n_min:n_max);

figure
hold on

plot(inverse_delta_t, y2_error_euler, '-o');
plot(inverse_delta_t, y2_error_midpoint, '-o');
plot(inverse_delta_t, y2_error_modified_euler, '-o')
plot(inverse_delta_t, y2_error_adams, '-o');
plot(inverse_delta_t, y2_error_runge_kutta, '-o');

hold off

set(gca, 'YScale', 'log', 'XScale', 'log')

xlabel('1/$\Delta$t')
ylabel('Absolute error')
title('Absolute Error at y(2) v. 1/$\Delta$t');
legend({'Euler', 'Midpoint', 'Modified Euler', '2-Step AB/AM', 'RK4'}, 'Location', 'bestoutside');

% saveas(gcf, 'results/q3-y2-errors.png')

% 
% Approximates y(t) given dy/dt and h.
% 
% Inputs:
%   dydt: derivative function y' = f(y, t)
%   y0: initial condition
%   y: exact solution function
%   a: lower t limit
%   b: upper t limit
%   n_min: largest step size, 2^(-n_min)
%   n_max: smallest step size, 2^(-n_max)
%   method_name: name of the approximation method

function y2_error = approximate(dydt, y0, y, a, b, n_min, n_max, method_name)
    switch method_name
        case 'euler'
            graph_title = "Euler's Method";
            method_func = @euler;
        case 'rk2'
            graph_title = "Midpoint Rule (Runge-Kutta 2nd Order)";
            method_func = @rk2;
        case 'eulerModified'
            graph_title = "Modified Euler's Method";
            method_func = @eulerModified;
        case 'adams2'
            graph_title = "2-Step A-B/A-M Predictor Corrector Scheme";
            method_func = @adams2;
        case 'rk4'
            graph_title = "Runge-Kutta 4th Order";
            method_func = @rk4;
        otherwise
            error('Unexpected method name.');
    end
    figure
    tiledlayout('vertical');

    % colors
    colormap spring
    c = spring(n_max - n_min + 1);
    colororder(c);

    ax_w = nexttile;
    ax_e = nexttile;

    title(ax_w, graph_title);
    subtitle(ax_w, 'y(t) v. t');
    ylabel(ax_w, 'y');
    xlabel(ax_w, 't');
    hold(ax_w, 'on');
    
    title(ax_e, 'Absolute Error v. t');
    ylabel(ax_e, 'Absolute error');
    xlabel(ax_e, 't');
    hold(ax_e, 'on');
    
    % absolute value of the error at y(2) versus 1/∆t on a log-log scale
    y2_error = zeros(n_max - n_min + 1, 1);

    for n = n_min : n_max
        % step size
        h = 2^(-n);
        t = a:h:b;

        % use the selected method to approximate y
        w = method_func(dydt, a, b, h, y0); 
        plot(ax_w, t, w);
        
        % plot errors
        e = abs(y(t).' - w);
        plot(ax_e, t, e);

        % save error of calculated w at t = 2
        y2_error(n - n_min + 1) = abs(w(2^(n+1) + 1) - y(2));
    end

    % plot the analytical solution in y(t) v. t
    t = a:1e-3:b;
    plot(ax_w, t, y(t), 'b')

    % colorbar 
    clim([n_min, n_max]);
    cb = colorbar(); 
    cb.Layout.Tile = 'south';
    cb.Label.String = 'n';
    cb.Ticks = n_min:1:n_max;
end




