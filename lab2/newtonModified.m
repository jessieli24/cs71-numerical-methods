%
% Jessie Li, CS 71 Fall 2023
% 
% Approximates root p of f using Modified Newton's method. 
% 
% Parameters:
% - f: function
% - df: f'
% - ddf: f''
% - x0: first initial approximation of root p
% - max_iter: maximum number of iterations
% - tol: tolerance
% 
% Returns:
% - n: number of iterations until dist(p_n, p_(n-1)) < tol
% - x: approximation of p after n iterations
% - err: array containing error after each iteration

function [n, x, err] = newtonModified(f, df, ddf, x0, max_iter, tol)
    err = zeros(max_iter, 1);
    
    for n = 1:max_iter
        x = x0 - (f(x0) * df(x0)) / ((df(x0))^2 - f(x0) * ddf(x0));
        
        err(n) = abs((x - x0)/x);
        if (err(n) < tol)
            break;
        end

        x0 = x;
    end

    fprintf('root: %0.5f, iterations: %d, error: %+10.6e\n', x, n, err(n));
end
