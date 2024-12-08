%
% Jessie Li, CS 71 Fall 2023
% 
% RK2 (midpoint rule) to approximate y(t) given dy/dt and h.
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

function w = rk2(f, a, b, h, y0)
     % number of subintervals
     n = (b-a) / h;

     w = zeros(size(y0, 1), n+1);
     w(:, 1) = y0;

     for i = 1 : n
         t = a + (i-1) * h;

         w_mid = w(:, i) + h/2 .* f(w(:, i), t);
         t_mid = t + h/2;
   
         w(:, i+1) = w(:, i) + h .* f(w_mid, t_mid);
     end
end