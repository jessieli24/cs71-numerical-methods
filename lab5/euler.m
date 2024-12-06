%
% Jessie Li, CS 71 Fall 2023
% 
% Euler's method to approximate y(t) given dy/dt and h.
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

function w = euler(dydt, a, b, h, y0)
    % number of subintervals
    n = (b-a) / h;
    
    w = zeros(n+1, 1);
    w(1) = y0;
    
    for i = 1 : n
        t = a + (i-1) * h;
        w(i+1) = w(i) + h * dydt(w(i), t);
    end
end
