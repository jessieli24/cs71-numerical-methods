% 
% main1b.m - solves the 1D model for heat transfer in tissue (from Lab 7) 
%            using Gauss-Seidel iteration. 
% 
% Jessie Li, CS 71 Fall 2023
%

function err = main1b()
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
    tolerances = [1e-4, 1e-5, 1e-6];
    maxIters = 2000;
    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    
    err = zeros(size(tolerances, 2), size(nSubintervals, 2));
    spectralRadEstimate = zeros(size(tolerances, 2), size(nSubintervals, 2));
    spectralRad = zeros(1, size(nSubintervals, 2));

    for t = 1 : size(tolerances, 2)
        tolerance = tolerances(t);
    
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
    
            % initial guess for Gauss-Seidel (analytic + random offset)
            rng(0)
            T0 = T_analytic(h:h:L-h).' + (-0.05 + 0.1*rand(n,1));
            
            % use Gauss-Seidel to solve the system for T = [T1 ... Tn]
            [T, spectralRadEstimate(t, j)] = gaussSeidel(A, b, T0, tolerance, maxIters);
    
            % plot T v. x
            xx = 0 : h : L;
            TT = [(Tc - Ta) T.' (Ts - Ta)];
            plot(xx, TT, '--', 'LineWidth', 2, 'DisplayName', sprintf('n + 1 = %d', n + 1));
    
            % calculate error
            err(t, j) = max(abs(TT - T_analytic(xx)));
    
            if t == 1
                % calculate the exact spectral radius
                spectralRad(j) = spectralRadius(A);
            end
      
        end
        hold off
    
        xlabel('x')
        ylabel('$\tilde{T}$')
        title({'Gauss-Seidel Steady-State Temperature', 'Distribution in the Absence of Heating'})
        subtitle(sprintf('tolerance = %1.0e', tolerances(t)))
        legend()
    end

    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    % plot ε versus (n + 1) on a log-log scale for various tolerances
    figure
    defaultColors()
    
    hold on
    for i = 1 : size(tolerances, 2)
        plot(nSubintervals, err(i, :), '-o', 'DisplayName', sprintf('tolerance = %1.0e', tolerances(i)));
    end
    hold off
    
    set(gca, 'YScale', 'log', 'XScale', 'log');
    ylabel('$\epsilon$')
    xlabel('n + 1')
    title({'Gauss-Seidel Error v. Number of Subintervals'})
    legend()

    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    % plot ρ versus (n + 1) for various tolerances
    figure
    defaultColors()
    
    hold on
    % exact
    plot(nSubintervals, spectralRad, '-o', 'DisplayName', 'Exact')
    
    % estimated
    for i = 1 : size(tolerances, 2)
        plot(nSubintervals, spectralRadEstimate(i, :), '--o', 'DisplayName', sprintf('tolerance = %1.0e', tolerances(i)));
    end
    hold off
    
    set(gca, 'YScale', 'log', 'XScale', 'log');
    ylabel('$\rho$')
    xlabel('n + 1')
    title({'Estimated Spectral Radius v. Number of Subintervals'})
    legend()   
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