%
% Jessie Li, CS 71 Fall 2023
% 
% RK4 to approximate y(t) given dy/dt and h.
% 
% Input:
%     f: derivative y' = f(y, t)
%     a: lower t limit
%     b: upper t limit
%     h: step size
%     y0: initial condition
% 
% Returns:
%     w: approximated y values

function w = rk4(f, a, b, h, y0)
    % number of subintervals
    n = (b-a) / h;

    w = zeros(size(y0, 1), n+1);
    w(:, 1) = y0;

    for i = 1 : n
        t = a + (i-1) * h;
        f1 = f(w(:, i), t);
        f2 = f(w(:, i) + h/2 .* f1, t + h/2);
        f3 = f(w(:, i) + h/2 .* f2, t + h/2);
        f4 = f(w(:, i) + h .* f3, t + h);
         
        w(:, i+1) = w(:, i) + h/6 .* (f1 + 2 .* (f2 + f3) + f4);
    end
end