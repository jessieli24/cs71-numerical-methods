%
% Jessie Li, CS 71 Fall 2023
% 
% 2-step Adams-Bashforth.
% 
% Input:
%     dydt: derivative y' = f(y, t)
%     a: lower t limit
%     b: upper t limit
%     h: step size
%     y0: initial condition
%     y: analytic solution as a function
% 
% Returns:
%     w: approximated y values

function w = ab2(dydt, a, b, h, y0, y)
    % number of subintervals
    n = (b-a) / h;

    % initial conditions
    w = zeros(n+1, 1);
    f = zeros(n+1, 1);

    w(1) = y0;
    f(1) = dydt(y0, a);
    
    % compute the first step using the analytic solution
    w(2) = y(a + h);

    for i = 2 : n
        t = a + (i-1) * h;
        f(i) = dydt(w(i), t);
  
        w(i+1) = w(i) + h/2 * (3 * f(i) - f(i-1));
    end

end