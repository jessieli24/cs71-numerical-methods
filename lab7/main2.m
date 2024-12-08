% 
% main2.m - solve a simple 1D model for heat transfer in tissue 
%           with temperature-dependent energy removal due to blood flow 
%           and superficial microwave heating
% 
% Jessie Li, CS 71 Fall 2023
%

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

% -------------------- constants ---------------------- %
% ----------------------------------------------------- %
L = 1;
lambda2 = 2.7;
Ta = 37;
Tc = 37;
Ts = 32;

T_analytic = @(x) (Ts - Tc) .* sinh(sqrt(lambda2) .* x) ./ sinh(sqrt(lambda2) * L);

n_subintervals = [5, 10, 20, 40, 80, 160];
errors = zeros(1, size(n_subintervals, 2));

% ----------------------------------------------------- %
% ----------------------------------------------------- %
figure 
defaultColors()

hold on

% plot analytic solution
xx = 0 : 1e-3 : L;
plot(xx, T_analytic(xx), 'magenta', 'LineWidth', 2, 'DisplayName', 'Analytic')

for j = 1 : size(n_subintervals, 2)
    n = n_subintervals(j) - 1;

    A = zeros(n, n);
    b = zeros(n, 1);
    h = L / (n + 1);

    for i = 1 : n
        A(i, i) = -(2 + h^2 * lambda2);
    
        if i < n
            A(i, i+1) = 1;
        end

        if i > 1
            A(i, i-1) = 1;
        end
    end
    
    % boundary conditions
    b(1) = -(Tc - Ta);
    b(n) = -(Ts - Ta);

    % solve the system for T = [T1 ... Tn]
    % T = A \ b;
    T = solveTridiagonal(A, b);

    % plot T v. x
    xx = 0 : h : L;
    TT = [(Tc - Ta) T.' (Ts - Ta)];

    plot(xx, TT, '--', 'LineWidth', 2, 'DisplayName', sprintf('n + 1 = %d', n + 1));

    % calculate error
    errors(j) = max(abs(TT - T_analytic(xx)));
    max(abs(TT - T_analytic(xx)))

end

hold off

xlabel('x')
ylabel('$\tilde{T}$')
title({'Steady-State Temperature Distribution', 'in the Absence of Heating'})
legend()

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% plot ε versus (n + 1) on a log-log scale
figure
defaultColors()

loglog(n_subintervals, errors, '-o')

ylabel('$\epsilon$')
xlabel('n + 1')
title({'Maximum Error v. Number of Subintervals', 'Without Microwave Heating'})

% centered difference approximations are O(h^2)
% --> halving subinterval should decrease error by a factor of 4

% ----------------------------------------------------- %
% ----------------------------------------------------- %
microwave_heating = 100;
gamma = 1 / L;

% small step size to approximate the exact solution
n_subintervals_exact = 16000;
n = n_subintervals_exact - 1;

A = zeros(n, n);
b = zeros(n, 1);
h = L / (n + 1);

xx_exact = 0 : h : L;

for i = 1 : n
    A(i, i) = -(2 + h^2 * lambda2);
    
    if i < n
        A(i, i+1) = 1;
    end

    if i > 1
        A(i, i-1) = 1;
    end

    b(i) = -h^2 * microwave_heating * exp(gamma * (L - xx_exact(i+1)));
end

% boundary conditions
b(1) = b(1) - (Tc - Ta);
b(n) = b(n) - (Ts - Ta);

% solve the system for T = [T1 ... Tn]
% T = A \ b;
T = solveTridiagonal(A, b);

% plot T v. x
TT_exact = [(Tc - Ta) T.' (Ts - Ta)];

plot(xx_exact, TT_exact, 'DisplayName', sprintf('n + 1 = %d', n + 1));

xlabel('x')
ylabel('$\tilde{T}$')
title({'"Exact" Steady-State Temperature Distribution', 'With Microwave Heating'})

% ----------------------------------------------------- %
% ----------------------------------------------------- %
errors_heat = zeros(1, size(n_subintervals, 2));

figure 
defaultColors()

hold on
plot(xx_exact, TT_exact, 'Color', [0.7 0 1], 'DisplayName', '"Exact"', 'LineWidth', 2);

for j = 1 : size(n_subintervals, 2)
    n = n_subintervals(j) - 1;

    A = zeros(n, n);
    b = zeros(n, 1);
    h = L / (n + 1);

    xx = 0 : h : L;

    for i = 1 : n
        A(i, i) = -(2 + h^2 * lambda2);
    
        if i < n
            A(i, i+1) = 1;
        end

        if i > 1
            A(i, i-1) = 1;
        end

        b(i) = -h^2 * microwave_heating * exp(gamma * (L - xx(i+1)));
    end
    
    % boundary conditions
    b(1) = b(1) - (Tc - Ta);
    b(n) = b(n) - (Ts - Ta);

    % solve the system for T = [T1 ... Tn]
    % T = A \ b;
    T = solveTridiagonal(A, b);

    % plot T v. x
    TT = [(Tc - Ta) T.' (Ts - Ta)];

    plot(xx, TT, 'DisplayName', sprintf('n + 1 = %d', n + 1));

    % calculate error
    errors_heat(j) = max(abs(TT - TT_exact(1 : n_subintervals_exact/n_subintervals(j) : end)));
    max(abs(TT - TT_exact(1 : n_subintervals_exact/n_subintervals(j) : end)))
end
hold off

xlabel('x')
ylabel('$\tilde{T}$')
title({'Steady-State Temperature Distribution', 'With Microwave Heating'})
legend()

% ----------------------------------------------------- %
% ----------------------------------------------------- %
% plot ε versus (n + 1) on a log-log scale
figure
defaultColors()

loglog(n_subintervals, errors_heat, '-o')

ylabel('$\epsilon$')
xlabel('n + 1')
title({'Maximum Error v. Number of Subintervals', 'with Microwave Heating'})

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
