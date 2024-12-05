% 
% Jessie Li, CS 71 Fall 2023
%
% Approximates root x of f using Newton's method. 
% Based on Algorithm 2.3 from Burden & Faires.
% 
% Parameters:
% - f: fixed-point function
% - df: derivative of f
% - x0: initial approximation of root
% - max_iter: maximum number of iterations
% - tol: tolerance
% 
% Returns:
% - n: number of iterations until dist(x_n, x_(n-1)) < tol
% - x: approximation of x after n iterations
% - err: array containing error after each iteration 

function [n, x, err] = newton(f, df, x0, max_iter, tol) 
    err = zeros(max_iter, 1);
    
    for n = 1:max_iter
        x = x0 - f(x0) / df(x0);
        err(n) = abs((x - x0) / x);
        
        if (err(n) < tol)
            break;
        end
    
        x0 = x;
    end

    fprintf('root: %0.5f, iterations: %d, error: %+10.6e\n', x, n, err(n));
end