%
% Jessie Li, CS 71 Fall 2023
% 
% Approximates root p of f using a variant of Newton's method that
% converges cubically.
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
% - n: number of iterations until dist(x_n, x_(n-1)) < tol
% - x: approximation of x after n iterations
% - err: array containing error after each iteration 

function [n, x, err] = newtonCubic(f, df, ddf, x0, max_iter, tol)
    err = zeros(max_iter, 1);

    for n = 1:max_iter
        q = f(x0);
        dq = df(x0);
        ddq = ddf(x0);

        x = x0 - q/dq - q^2 * ddq/(2*dq^3);
        
        err(n) = abs((x - x0)/x);
        if (err(n) < tol)
            break;
        end

        x0 = x;
    end

    fprintf('root: %0.5f, iterations: %d, error: %+10.6e\n', x, n, err(n));
end
