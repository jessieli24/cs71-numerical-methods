%
% Jessie Li, CS 71 Fall 2023
% 
% 4-step Adams-Bashforth/Adams-Moulton predictor corrector scheme
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

function w = adams4(dydt, a, b, h, y0)
    % number of subintervals
    n = (b-a) / h;

    w = zeros(size(y0, 1), n+1);
    f = zeros(size(y0, 1), n+1);

    % initial conditions
    w(:, 1) = y0;
    f(:, 1) = dydt(y0, a);
    
    % compute the first three (next) steps using RK4
    for i = 1:3
        t = a + (i-1) * h;
        
        f1 = dydt(y0, t);
        f2 = dydt(y0 + h/2 .* f1, t + h/2);
        f3 = dydt(y0 + h/2 .* f2, t + h/2);
        f4 = dydt(y0 + h .* f3, t + h);

        w(:, i+1) = y0 + h/6 .* (f1 + 2 .* (f2 + f3) + f4);
        f(:, i+1) = dydt(w(:, i+1), t+h);

        y0 = w(:, i+1);
    end

    for i = 4 : n
        t = a + i * h;
        
        % predict w(i+1) with the Adams-Bashforth 4-step predictor
        w(:, i+1) = w(:, i) + h/24 .* (55 .* f(:, i)- 59 .* f(:, i-1) + 37 .* f(:, i-2) - 9 .* f(:, i-3));
        f(:, i+1) = dydt(w(:, i+1), t);
        
        % correct w(i+1) with the Adams-Moulton 4-step corrector
        w(:, i+1) = w(:, i) + h/720 .* (251 .* f(:, i+1) + 646 .* f(:, i) - 264 .* f(:, i-1) + 106 .* f(:, i-2) - 19 .* f(:, i-3));
        f(:, i+1) = dydt(w(:, i+1), t);
    end
end
