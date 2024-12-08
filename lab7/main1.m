% 
% main1.m - solves a BVP using a nonlinear shooting method with a 4-step
%           Adams predictor-corrector for y(x), the deflection of a 
%           uniformly loaded, long rectangular beam under an axial 
%           tension force
% 
% Jessie Li, CS 71 Fall 2023
%

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

% -------------------- constants ---------------------- %
% ----------------------------------------------------- %
L = 50;
D = 8.5 * 1e7;
S = 100;
q = 1000;

% decompose y'' = f(x, y, y') into a linear system
F = @(z, x) [z(2); 
            (S/D * z(2) + q*x / (2*D) * (x-L) * z(1)) * (1 + z(2)^2) ^ 1.5];

% endpoints, x = a and x = b
a = 0;
b = 50;

% boundary conditions, y(a) and y(b)
ya = 0;
yb = 0;

% parameters for Newton's method
max_iterations = 20;
tolerance = 1e-14;

% parameters for 4-step Adams
n = 5000;
h = (b-a)/n; % = 0.01

% initial guess for y'(a)
u = -1.5;

% δf/δy 
df_dy = @(w, x) q*x / (2*D) * (x-L) * (1 + w(2)^2) ^ 1.5;

% δf/δy'
df_dydt = @(w, x) 3 * (1 + w(2)^2) ^ (1/2) * w(2) ...
                    * (S/D * w(2) + q*x / (2*D) * (x - L) * w(1)) ...
                    + ((1 + w(2)^2) ^ 1.5) * S/D;

% decompose g'' = δf/δy * g + δf/δy' * g'
G = @(v, w, x) [v(2); 
                df_dy(w, x) * v(1) +  df_dydt(w, x) * v(2)];

% ----------------------------------------------------- %
% solve BVP
% ----------------------------------------------------- %
% results
result_message = '';
newton_errors = zeros(1, max_iterations);

xx = a : h : b;

figure
gradientBlueColors(9)
hold on

yline(0, '--')

for k = 1 : max_iterations 
    % solve y'' = f(x, y, y') with y(a) = α, y'(a) = u0 for z1(b, u_n)
    w = zeros(2, n + 1);      
    f = zeros(2, n + 1);

    % solve g'' = δf/δy * g + δf/δy' * g' with g(a) = 0, g'(a) = 1 for g(b, u_n)
    v = zeros(2, n + 1);
    g = zeros(2, n + 1);

    % solve system to get z1(b, u) then δz1(b, u)/δu using 4-step Adams
    w(:, 1) = [ya; u];
    f(:, 1) = F(w(:, 1), a);

    v(:, 1) = [0; 1];
    g(:, 1) = G(v(:, 1), w(:, 1), a);

    % compute the first three steps using RK4
    for i = 1 : 3
        t = a + (i-1) * h;

        % RK4 for [z1(x, u_n); z2(x, u_n)]
        f1 = F(w(:, i), t);
        f2 = F(w(:, i) + h/2 .* f1, t + h/2);
        f3 = F(w(:, i) + h/2 .* f2, t + h/2);
        f4 = F(w(:, i) + h .* f3, t + h);

        w(:, i+1) = w(:, i) + h/6 .* (f1 + 2 .* (f2 + f3) + f4);
        f(:, i+1) = F(w(:, i+1), t + h);

        % RK4 for [v1(x, u_n); v2(x, u_n)]
        g1 = G(v(:, i), w(:, i), t);
        g2 = G(v(:, i) + h/2 .* g1, w(:, i), t + h/2);
        g3 = G(v(:, i) + h/2 .* g2, w(:, i), t + h/2);
        g4 = G(v(:, i) + h .* g3, w(:, i), t + h);

        v(:, i+1) = v(:, i) + h/6 .* (g1 + 2 .* (g2 + g3) + g4);
        g(:, i+1) = G(v(:, i+1), w(:, i+1), t + h);
    end
        
    % 4-step Adams predictor-corrector
    for i = 4 : n 
        t = a + i * h;

        % predict w(i+1) with the Adams-Bashforth 4-step predictor
        w(:, i+1) = w(:, i) + h/24 .* (55 .* f(:, i) - 59 .* f(:, i-1) + 37 .* f(:, i-2) - 9 .* f(:, i-3));
        f(:, i+1) = F(w(:, i+1), t);

        v(:, i+1) = v(:, i) + h/24 .* (55 .* g(:, i) - 59 .* g(:, i-1) + 37 .* g(:, i-2) - 9 .* g(:, i-3));
        g(:, i+1) = G(v(:, i+1), w(:, i+1), t);

        % correct w(i+1) with the Adams-Moulton 4-step corrector
        w(:, i+1) = w(:, i) + h/720 .* (251 .* f(:, i+1) + 646 .* f(:, i) - 264 .* f(:, i-1) + 106 .* f(:, i-2) - 19 .* f(:, i-3));
        f(:, i+1) = F(w(:, i+1), t);

        v(:, i+1) = v(:, i) + h/720 .* (251 .* g(:, i+1) + 646 .* g(:, i) - 264 .* g(:, i-1) + 106 .* g(:, i-2) - 19 .* g(:, i-3));
        g(:, i+1) = G(v(:, i+1), w(:, i+1), t);     
    end

    newton_errors(k) = abs(w(1, end) - yb);

    % stop if error is less than tolerance
    if abs(w(1, end) - yb) < tolerance
        % plot final solution for y(x)
        p = plot(xx, w(1, :), 'magenta');
        scatter(50, w(1, end), 'filled', 'MarkerFaceColor', 'magenta');
        result_message = sprintf('found u = %.5f after %d iterations', u, k);
        break;
    
    elseif k == max_iterations
        % plot y(x)
        p = plot(xx, w(1, :));
        scatter(50, w(1, end), 'filled', 'MarkerFaceColor', get(p, 'Color'));
        result_message = sprintf('exceeded maximum number of iterations, error: %f', abs(w(1, end) - yb));
        break;
    
    end

    % plot y(x)
    p = plot(xx, w(1, :));
    scatter(50, w(1, end), 'filled', 'MarkerFaceColor', get(p, 'Color'));

    % update u with Newton's method
    u = u - (w(1, end) - yb) / v(1, end);
end

% boundary conditions, y(a) and y(b)
scatter([0, 50], [0, 0], 'filled', 'MarkerFaceColor', 'magenta');
hold off

xlabel('x')
ylabel('y')
title({'Deflection of a Uniformly Loaded Beam', 'Under an Axial Tension Force'})

fprintf(result_message);

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% final solution for y as a function of x
xx = 0 : h : 50;

figure
defaultColors()

hold on
plot(xx, w(1, :));  % y(x)

% boundary conditions, y(a) and y(b)
scatter([0, 50], [0, 0], 'filled');
hold off

xlabel('x')
ylabel('y')

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% final solution for y' as a function of x
figure
defaultColors()

plot(xx, w(2, :));  % y'(x)

xlabel('x')
ylabel('y''')

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% convergence of Newton's method
figure
defaultColors()

% plot error at x = b v. iteration (Newton's)
semilogy(1 : k, newton_errors(1 : k), '-o')

xlabel('Iteration')
ylabel('Absolute error of y at x = b')
title('Convergence of Newton''s Method')

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% order of convergence of Newton's method
alpha = zeros(1, k - 2);

for i = 2 : k-1
    alpha(i-1) = (log(newton_errors(i+1)) - log(newton_errors(i))) / ...
                 (log(newton_errors(i)) - log(newton_errors(i-1)));
end

figure
plot(2 : k-1, alpha, '-o')

xlabel('Iteration')
ylabel('Order of convergence $\alpha$')
title('Order of Convergence of Newton''s Method')

% ----------------------------------------------------- %
% find y(x) and y'(x) using 4-Step Adams with various step sizes with 
% the calculated value of u
% ----------------------------------------------------- %
a = 0;
b = 50;

stepSizes = [10, 1, 0.1, 0.01, 0.001, 0.0001];
w = zeros(size(stepSizes, 2) * 2, (b-a)/min(stepSizes));

for i = 1 : 2 : size(stepSizes, 2) * 2 - 1
    h = stepSizes((i+1)/2);
    n = (b-a)/h;
    w(i:i+1, 1:n+1) = adams4(F, a, b, h, [ya; u]);
end

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% plot y(x) for each step size
figure
gradientCoolColors(size(stepSizes, 2))

hold on
yline(0, '--', 'DisplayName', 'y = 0')

for i = 1 : 2 : size(stepSizes, 2) * 2 - 1
    h = stepSizes((i+1)/2);
    xx = a : h : b;
    n = (b-a)/h;
    
    p = plot(xx, w(i, 1:n+1), 'DisplayName', sprintf('h = %.4g', h));
end

% boundary conditions, y(a) and y(b)
scatter([0, 50], [0, 0], 'filled', 'MarkerFaceColor', get(p, 'Color'), 'DisplayName', 'B.C.');
hold off

xlabel('x')
ylabel('y')
title({'4-Step Adams Predictor-Corrector', 'with Various Step Sizes'})
legend("FontSize", 14, "Location", "southeast")

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% plot y'(x) for each step size
figure
gradientCoolColors(size(stepSizes, 2))

hold on
for i = 2 : 2 : size(stepSizes, 2) * 2
    h = stepSizes(i/2);
    xx = a : h : b;
    n = (b-a)/h;
    
    plot(xx, w(i, 1 : n+1), 'DisplayName', sprintf('h = %.4g', h))
end

xlabel('x')
ylabel('y''')
title({'4-Step Adams Predictor-Corrector', 'with Various Step Sizes'})
legend('Location', 'southeast', 'FontSize', 14)

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% for each step size, plot the truncation error in y of the 4-step 
% Adams predictor-corrector error v. x assuming smallest step size 
% is exact solution

figure
gradientBlueColors(size(stepSizes, 2))

hold on
for i = 1 : 2 : size(stepSizes, 2) * 2 - 2
    h = stepSizes((i+1)/2);
    xx = a : h : b;
    n = (b-a)/h;

    error = abs(w(i, :) - w(end-1, :));
    plot(xx, error(1 : n+1), 'DisplayName', sprintf('h = %.4g', h))
end
hold off

set(gca, 'YScale', 'log')
xlabel('x')
ylabel('Absolute error in y')
title({'Absolute Error in y v. x for', 'Various Step Sizes'})
legend()

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% error in y' v. x for each step size
figure
gradientBlueColors(size(stepSizes, 2))

hold on
for i = 2 : 2 : size(stepSizes, 2) * 2 - 1
    h = stepSizes(i/2);
    xx = a : h : b;
    n = (b-a)/h;

    error = abs(w(i, :) - w(end, :));
    plot(xx, error(1 : n+1), 'DisplayName', sprintf('h = %.4g', h))
end
hold off

set(gca, 'YScale', 'log')
xlabel('x')
ylabel('Absolute error in y''')
title({'Absolute Error in y'' v. x for', 'Various Step Sizes'})
legend()

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% plot error of y at x = b v. step size
figure
defaultColors()

errors_yb = zeros(size(stepSizes));

for i = 1 : 2 : size(stepSizes, 2) * 2 - 2
    h = stepSizes((i+1)/2);
    n = (b-a)/h;
    
    % assume smallest step size is "exact"
    errors_yb((i+1)/2) = abs(w(end-1, end) - w(i, n+1));
end

loglog(1 ./ stepSizes, errors_yb, '-o')

xlabel('1/h')
ylabel('Absolute error of y(b)')
title({'Absolute Error of y at x = b versus 1/h'})

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% plot error of y' at x = b v. step size
figure
defaultColors()

errors_dyb = zeros(1, size(stepSizes, 2) - 1);

for i = 2 : 2 : size(stepSizes, 2) * 2 - 2
    h = stepSizes(i/2);
    n = (b-a)/h;

    % assume smallest step size is "exact"
    errors_dyb(i/2) = abs(w(end, end) - w(i, n+1));
end

loglog(1 ./ stepSizes(1:end-1), errors_dyb, '-o')

set(gca, 'YScale', 'log')
xlabel('1/h')
ylabel('Absolute error of y''(b)')
title({'Absolute Error of y'' at x = b versus 1/h'})

% ---------------- helper functions ------------------- %
% ----------------------------------------------------- %

function defaultColors()
    color_order = [0.37 0.60 0.94
                   0.05 0.26 0.57
                   0.98 0.58 0.89
                   0.99 0.82 0.54
                   0.81 0.59 0.95
                   0.53 0.98 0.84];
    
    colororder(color_order)
end

function gradientBlueColors(i)
    colormap sky
    c = sky(i);
    colororder(c);
end

function gradientCoolColors(i)
    colormap cool
    c = cool(i);
    colororder(c);
end


