% 
% main1c.m - implements the Crank-Nicholson time-stepping algorithm to 
%           account for the time derivative term in the bioheat equation
% 
% Jessie Li, CS 71 Fall 2023
%

function err = main1c()
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
    tolerance = 1e-8;
    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    err = zeros(size(nSubintervals));
    
    figure 
    defaultColors()
    
    hold on
    % plot analytic solution
    xx = 0 : 1e-3 : L;
    plot(xx, T_analytic(xx), 'magenta', 'LineWidth', 2, 'DisplayName', 'Analytic')
    
    for j = 1 : size(nSubintervals, 2)
        n = nSubintervals(j) - 1;
    
        dt = 1/(4*n);        % ∆t/h << 1/2
        maxIters = 2/dt;     % time = 2
        
        [A, B, c, h] = getCrankNicholsonMatrix(n, dt);
        
        A = LUDecomposeTridiagonal(A);
        T0 = zeros(n, 1) - 5;
    
        for k = 1 : maxIters
            % solve the system A * T(k+1) = B * T(k) + c
            b = B * T0 + c;
            T1 = forwardSubstitute(A, b);
            T1 = backwardSubstitute(A, T1);
            
            % check if steady-state
            if max(abs(T1 - T0)) < tolerance
                fprintf(sprintf('reached steady-state after %d iterations\n', k));
                break;
            elseif k == maxIters
                fprintf(sprintf('reached maximum number of iterations: %d\n', k));
            end
            
            T0 = T1;
        end
        
        % plot the steady-state solution, T v. x
        xx = 0 : h : L;
        TT = [(Tc - Ta) T1.' (Ts - Ta)];
    
        plot(xx, TT, '--', 'LineWidth', 2, 'DisplayName', sprintf('n + 1 = %d', n + 1));
    
        % calculate error
        err(j) = max(abs(TT - T_analytic(xx)));
    end
    hold off
    
    xlabel('x')
    ylabel('$\tilde{T}$')
    title({'Crank-Nicholson Steady-State Temperature', 'Distribution in the Absence of Heating'})
    legend()

    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    % plot ε versus (n + 1) on a log-log scale
    figure
    defaultColors()
    
    loglog(nSubintervals, err, '-o')
    
    ylabel('$\epsilon$')
    xlabel('n + 1')
    title({'Crank-Nicholson Error v. Number of Subintervals'})
end

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