% 
% main1a.m - solves the 1D model for heat transfer in tissue (from Lab 7) 
%            with Crout LU decomposition. 
% 
% Jessie Li, CS 71 Fall 2023
%

function err = main1a()
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
    
    nSubintervals = [5, 10, 20, 40, 80, 160];
    
    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    err = zeros(1, size(nSubintervals, 2));
    
    figure 
    defaultColors()
    
    hold on
    % plot analytic solution
    xx = 0 : 1e-3 : L;
    plot(xx, T_analytic(xx), 'magenta', 'LineWidth', 2, 'DisplayName', 'Analytic')
    
    for j = 1 : size(nSubintervals, 2)
        n = nSubintervals(j) - 1;
    
        % construct system of finite difference equations
        [A, b, h] = getFiniteDiffMatrix(n);
    
        % solve the system for T = [T1 ... Tn]
        T = solveTridiagonal(A, b);
    
        % plot T v. x
        xx = 0 : h : L;
        TT = [(Tc - Ta) T.' (Ts - Ta)];
    
        plot(xx, TT, '--', 'LineWidth', 2, 'DisplayName', sprintf('n + 1 = %d', n + 1));
    
        % calculate error
        err(j) = max(abs(TT - T_analytic(xx)));
        max(abs(TT - T_analytic(xx)))
    
    end
    hold off
    
    xlabel('x')
    ylabel('$\tilde{T}$')
    title({'Crout LU Steady-State Temperature', 'Distribution in the Absence of Heating'})
    legend()
    
    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    % plot ε versus (n + 1) on a log-log scale
    figure
    defaultColors()
    
    loglog(nSubintervals, err, '-o')
    
    ylabel('$\epsilon$')
    xlabel('n + 1')
    title({'Crout LU Error v. Number of Subintervals'})
end

% ----------------------------------------------------- %
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
