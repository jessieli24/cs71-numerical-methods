% 
% main2b.m - numerical differentiation
% 
% Jessie Li, CS 71 Fall 2023
%

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

% -------------------- constants ---------------------- %
% ----------------------------------------------------- %
TOL = 1e-10;
MAX_ITER = 20;

CG = 1.25;
GF = 1.26;
EF = 1.82;
CE = 2.41;

w_rad_min = 550;      % radians per minute
w_rad_sec = 550/60;   % radians per second
% ----------------- solve for phi --------------------- %
% ----------------------------------------------------- %
phi = main2a(true);

% ----------------- solve for beta -------------------- %
% ----------------------------------------------------- %

% store [β θ3]
beta = zeros(362, 1);
t3 = zeros(362, 1);

% initial guess for θ = 0
beta(1) = deg2rad(78);
t3(1) = deg2rad(218.5);

for t = 1 : 361
    t4 = phi(t) + deg2rad(149 + 180);

    f = @(x) [GF * cos(x(1)) + EF * cos(x(2)) + CE * cos(t4) - CG;
              GF * sin(x(1)) + EF * sin(x(2)) + CE * sin(t4)];
    
    J = @(x) [-GF*sin(x(1)), -EF*sin(x(2));
              GF*cos(x(1)),  EF*cos(x(2))];
    
    % use previously found solution as initial starting guess
    x0 = [beta(t); t3(t)];
    [~, x, ~] = newtonNonlinear(f, J, x0, MAX_ITER, TOL);
    beta(t+1) = mod(x(1), 2*pi);
    t3(t+1) = x(2);
end

beta = beta(2:end);

% plot β v. θ
f = figure;
plot(deg2rad(0:1:360), beta);

xTicksRadians;
yTicksRadians;

title('$\beta$ v. $\theta$');
xlabel('$\theta$');
ylabel('$\beta$');

% ---------------- first derivative ------------------- %
% ----------------------------------------------------- %

% dβ/dθ with a forward difference 
dbeta_dtheta_forward = forwardDifference1(beta, deg2rad(1));

% dβ/dθ with a centered difference
dbeta_dtheta_centered = centeredDifference1(beta, deg2rad(1));

% plot both curves on the same graph
plotDerivatives(deg2rad(0:1:359), dbeta_dtheta_forward, ...
                deg2rad(1:1:359), dbeta_dtheta_centered, ...
                '$\theta$', '$d\phi/d\theta$');

% plot difference between forward and centered
plotDifference(dbeta_dtheta_forward(2:end), ...
               dbeta_dtheta_centered, ...
               deg2rad(1:1:359), '$\theta$', '$d\beta/d\theta$');

% ---------------- second derivative ------------------ %
% ----------------------------------------------------- %

% d2β/dθ2 with a second forward difference 
d2beta_dtheta2_forward = forwardDifference2(beta, deg2rad(1));

% d2β/dθ2 with a centered difference
d2beta_dtheta2_centered = centeredDifference2(beta, deg2rad(1));

% plot d2β/dθ2 curves on the same graph
plotDerivatives(deg2rad(0:1:358), d2beta_dtheta2_forward, ...
                deg2rad(1:1:359), d2beta_dtheta2_centered, ...
                '$\theta$', '$d^2\beta/d\theta^2$');

% plot differences between forward and centered for d2β/dθ2
plotDifference(d2beta_dtheta2_forward(2:end), ...
               d2beta_dtheta2_centered(1:end-1), ...
               deg2rad(1:1:358), '$\theta$', '$d^2\beta/d\theta^2$');

% ------------------- in rad/sec ---------------------- %
% ----------------------------------------------------- %

t = deg2rad(0:1:360) / w_rad_sec;

% plot β v. t
f = figure;
plot(t, beta);

yTicksRadians;
title('$\beta$ v. $t$');
xlabel('$t$');
ylabel('$\beta$');

% forward dβ/dt
dbeta_dt_forward = dbeta_dtheta_forward * w_rad_sec;

% centered dβ/dt
dbeta_dt_centered = dbeta_dtheta_centered * w_rad_sec;

% plot dβ/dt curves on the same graph
plotDerivatives(t(1:360), dbeta_dt_forward, ...
                t(2:360), dbeta_dt_centered, ...
                '$t$', '$d\beta/dt$');
                  
% plot difference between forward and centered for dβ/dt
plotDifference(dbeta_dtheta_forward(2:end), ...
               dbeta_dtheta_centered, ...
               t(2:360), '$t$', '$d\beta/dt$');

% forward d2β/dt2
d2beta_dt2_forward = d2beta_dtheta2_forward * w_rad_sec^2;

% centered d2β/dt2
d2beta_dt2_centered = d2beta_dtheta2_centered * w_rad_sec^2;

% plot d2β/dt2 curves on the same graph
plotDerivatives(t(1:359), d2beta_dt2_forward, ...
                      t(2:360), d2beta_dt2_centered, ...
                      '$t$', '$d^2\beta/dt^2$');
                  
% plot differences between forward and centered for d2β/dt2
plotDifference(d2beta_dt2_forward(2:end), ...
                d2beta_dt2_centered(1:end-1), ...
                t(2:359), '$t$', '$d^2\beta/dt^2$');

% ---------------- helper functions ------------------- %
% ----------------------------------------------------- %
function f = plotDerivatives(xf, yf, xc, yc, xLabelText, yLabelText)
    f = figure;
    plot(xf, yf, xc, yc);
    
    xlabel(xLabelText);
    ylabel(yLabelText);
    legend({'forward', 'centered'});
    xTicksRadians;
end

function f = plotDifference(forward, centered, tt, xLabelText, yLabelText)
    delta = forward - centered;
    
    deltaPositive = delta;
    deltaPositive(delta < 0) = NaN;
    deltaNegative = delta;
    deltaNegative(delta >= 0) = NaN;
    
    f = figure;
    plot(tt, deltaPositive, tt, abs(deltaNegative), 'r');
    set(gca, 'YScale', 'log');
    
    xlabel(xLabelText);
    ylabel('Absolute difference');
    title({'Absolute difference between forward and', sprintf('centered estimates of %s', yLabelText)});
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
