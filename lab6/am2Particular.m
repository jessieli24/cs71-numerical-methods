%
% Jessie Li, CS 71 Fall 2023
% 
% 2-step Adams-Moulton, specific to y' = -7y (in Q2)
% 
% Input:
%     a: lower t limit
%     b: upper t limit
%     h: step size
%     y0: initial condition
%     y: analytic solution as a function
% 
% Returns:
%     w: approximated y values

function w = am2Particular(a, b, h, y0, y)
    % number of subintervals
    n = (b-a) / h;

    w = zeros(n+1, 1);
    w(1) = y0;
    
    % compute the first step using the analytic solution
    w(2) = y(a + h);

    for i = 2 : n
        % *note: specific to y' = -7y
        w(i+1) = 1/(1 + 35/12 * h) * ((1 - 56/12 * h) * w(i) + 7/12 * h * w(i-1));
    end
end
