%
% Jessie Li, CS 71 Fall 2023
% 
% Approximates root x of f using the Secant method. 
% Based on Algorithm 2.4 from Burden & Faires.
% 
% Parameters:
% - f: function
% - x0: first initial approximation of root x
% - x1: second initial approximation of root x
% - max_iter: maximum number of iterations
% - tol: tolerance
% 
% Returns:
% - n: number of iterations until dist(x_n, x_(n-1)) < tol
% - x: approximation of x after n iterations
% - err: array containing error after each iteration

function [n, x, err] = secant(f, x0, x1, max_iter, tol)
    err = zeros(max_iter, 1);
    q0 = f(x0);
    q1 = f(x1);
    
    for n = 1:max_iter
        x = x1 - q1 * (x1 - x0) / (q1 - q0);
        
        err(n) = abs((x - x1)/x);
        if (err(n) < tol)
            break;
        end
    
        x0 = x1;
        q0 = q1;
        x1 = x;
        q1 = f(x);
    end
    
    fprintf('root: %0.5f, iterations: %d, error: %+10.6e', x, n, err(n));
end
