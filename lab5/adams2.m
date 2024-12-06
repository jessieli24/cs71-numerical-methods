%
% Jessie Li, CS 71 Fall 2023
% 
% 2-step Adams-Bashforth/Adams-Moulton predictor corrector scheme
% with a single correction.
% 
% Input:
%     dydt: derivative y' = f(y, t)
%     a: lower t limit
%     b: upper t limit
%     h: step size
%     y0: initial condition
% 
% Returns:
%     w: approximated y values

function w = adams2(dydt, a, b, h, y0)
    % number of subintervals
    n = (b - a) / h;

    w = zeros(n+1, 1);
    f = zeros(n+1, 1);

    % initial conditions
    w(1) = y0;
    f(1) = dydt(y0, a);
    
    % compute the first step using RK-4
    f1 = dydt(y0, a);
    f2 = dydt(y0 + h/2 * f1, a + h/2);
    f3 = dydt(y0 + h/2 * f2, a + h/2);
    f4 = dydt(y0 + h * f3, a + h);

    w(2) = y0 + h/6 * (f1 + 2 * (f2 + f3) + f4);
    
    for i = 2 : n
        t = a + (i-1) * h;
        f(i) = dydt(w(i), t);
        
        % predict w(i+1) with the Adams-Bashforth 2-step predictor
        w(i+1) = w(i) + h/2 * (3 * f(i) - f(i-1));
        f(i+1) = dydt(w(i+1), t + h);
        
        % correct w(i+1) with the Adams-Moulton 2-step corrector
        w(i+1) = w(i) + h/12 * (5 * f(i+1) + 8 * f(i) - f(i-1));
    end
end