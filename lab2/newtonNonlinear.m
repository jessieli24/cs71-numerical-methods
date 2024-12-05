
%
% Jessie Li, CS 71 Fall 2023
% 
% Approximates root x of f using Newton's method for nonlinear systems. 
% 
% Parameters:
% - f: nx1 array of functions
% - J: Jacobian of f
% - x0: initial approximation of x
% - max_iter: maximum number of iterations
% - tol: tolerance
% 
% Returns:
% - n: number of iterations until dist(x_n, x_(n-1)) < tol
% - x: approximation of x after n iterations
% - err: array containing error after each iteration

function [n, x, err] = newtonNonlinear(f, J, x0, max_iter, tol)
    err = zeros(max_iter, 1);
    
    for n = 1:max_iter
        y = J(x0)\f(x0);   % solves Jy = f
        
        x = x0 - y;
        
        err(n) = norm(x - x0)/norm(x);
        if (err(n) < tol)
            break;
        end

        x0 = x;
    end

    fprintf('root: %s, iterations: %d, error: %+10.6e\n', mat2str(x, 5), n, err(n));
    fprintf('root, in degrees: %s\n', mat2str(rad2deg(x)));
end