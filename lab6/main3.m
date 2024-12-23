% 
% main3.m - solves the Van der Pol equation 
%           y'' = a*y' - (y')^3 - y
%
% Jessie Li, CS 71 Fall 2023
% 
% Inputs:
%   methodName (string): name of the method to use
%

function main3(methodName)
    % -------------------- constants ---------------------- %
    % ----------------------------------------------------- %
    tMin = 0;
    tMax = 100;
    
    aMin = 0.5;
    aMax = 5.5;
    % 0.5, 1.5, 2.5, 3.5, 4.5 and 5.5
    aRange = aMin : 1 : aMax;
    
    % y(0) = 0.0, y'(0) = 0.1
    z0 = [0.0; 0.1];

    switch methodName
        case 'rk2'
            graphTitle = 'Runge-Kutta 2nd Order';
            methodFunc = @rk2;
        case 'adams4'
            graphTitle = '4-Step Adams Predictor-Corrector';
            methodFunc = @adams4;
        otherwise
            warning('Method not implemented, using 4-Step Adams Predictor-Corrector');
            graphTitle = '4-Step Adams Predictor-Corrector';
            methodFunc = @adams4;
    end
    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    % preview the solution with rk2
    figure
    colormap spring
    c = spring(length(aRange));
    colororder(c);
    
    hold on
    for a = aRange
        
        % decompose into a linear system: z1 = y, z2 = y'
        f = @(z, t) [z(2); 
                     a * z(2) - z(2)^3 - z(1)];
    
        w = methodFunc(f, tMin, tMax, 1e-3, z0);
    
        plot(w(1, :), w(2, :))
    end
    hold off
    
    xlabel('Voltage (y)')
    ylabel('Current (y'')')
    title(graphTitle)
    
    % colorbar 
    clim([aMin, aMax]);
    cb = colorbar(); 
    cb.Label.String = 'a';
    cb.Ticks = aRange;
end
