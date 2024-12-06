% 
% main2a.m - numerical differentiation
% 
% Jessie Li, CS 71 Fall 2023
%
% Inputs:
%   phiOnly (boolean): if true, stop after calculating phi

function phi = main2a(phiOnly)
    set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'DefaultTextInterpreter', 'latex');
    set(groot, 'DefaultLegendInterpreter', 'latex');
    
    % -------------------- constants ---------------------- %
    % ----------------------------------------------------- %
    TOL = 1e-10;
    MAX_ITER = 20;
    
    DC = 7;
    CB = 2.36;
    AB = 6.86;
    DA = 1.84;
    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    
    % store [φ θ]
    phi = zeros(362, 1);
    t3 = zeros(362, 1);
    
    % initial guess for θ = 0
    phi(1) = deg2rad(28.5);
    t3(1) = deg2rad(208);
    
    for t = 0 : 360
        t4 = deg2rad(t + 180);
    
        f = @(x) [CB * cos(x(1)) + AB * cos(x(2)) + DA * cos(t4) - DC;
                  CB * sin(x(1)) + AB * sin(x(2)) + DA * sin(t4)];
        
        J = @(x) [-CB*sin(x(1)), -AB*sin(x(2));
                  CB*cos(x(1)),  AB*cos(x(2))];
        
        % use previously found solution as initial starting guess
        x0 = [phi(t+1); t3(t+1)];
        [~, x, ~] = newtonNonlinear(f, J, x0, MAX_ITER, TOL);
        phi(t+2) = mod(x(1), 2*pi);
        t3(t+2) = x(2);
    end
    
    phi = phi(2:end);
    if (phiOnly)
        return;
    end

    % plot φ v. θ
    f = figure;
    plot(deg2rad(0:1:360), phi);
    
    xTicksRadians;
    yTicksRadians;
    
    title('$\phi$ v. $\theta$');
    xlabel('$\theta$');
    ylabel('$\phi$');
    
    % dφ/dθ with a forward difference 
    dphi_dtheta_forward = forwardDifference1(phi, deg2rad(1));
    
    % dφ/dθ with a centered difference
    dphi_dtheta_centered = centeredDifference1(phi, deg2rad(1));
    
    % plot both curves on the same graph
    f = figure;
    plot(deg2rad(0:1:359), dphi_dtheta_forward, ...
         deg2rad(1:1:359), dphi_dtheta_centered);
    
    xlabel('$\theta$');
    ylabel('$d\phi/d\theta$');
    legend({'forward', 'centered'});
    xTicksRadians;
    
    % plot difference between forward and centered
    delta = dphi_dtheta_forward(2:end) - dphi_dtheta_centered;
    
    deltaPositive = delta;
    deltaPositive(delta < 0) = NaN;
    deltaNegative = delta;
    deltaNegative(delta >= 0) = NaN;
    
    f = figure;
    xx = deg2rad(1:1:359);
    plot(xx, deltaPositive, xx, abs(deltaNegative), 'r');
    set(gca, 'YScale', 'log');
    
    xlabel('$\theta$');
    ylabel('Absolute difference');
    title('Absolute difference between forward and centered estimates of $d\phi/d\theta$');
    legend({'forward $\geq$ centered', 'forward $<$ centered'});
    xTicksRadians;
end

function xTicksRadians
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2' '2\pi'})
end

function yTicksRadians
    yticks([0 pi/4 pi/2 3*pi/4 pi, 5*pi/4, 3*pi/2, 2*pi])
    yticklabels({'0','\pi/4','\pi/2','3\pi/4' '\pi', '5\pi/4', '3\pi/2', '2\pi'})
end