%
% Jessie Li, CS 71 Fall 2023
% 
% Approximates a definite integral with Romberg's method.
% 
% Input:
%     f: function
%     a: lower limit
%     b: upper limit
%     tol: tolerance
%     max_steps: maximum number of iterations
% 
% Returns:
%     y: Romberg table
%     n: number of times subintervals were halved
%     evals: number of function evaluations

function [y, n, evals] = romberg(f, a, b, tol, max_steps)
    y = zeros(max_steps+1, max_steps+1);
    evals = 0;

    for n = 0 : max_steps
        for k = 0 : n
       
            % initialize with a composite trapezoidal approximation
            if k == 0
                [y(n+1, k+1), trapEvals] = compositeTrapezoidal(f, a, b, 2^n);
                evals = evals + trapEvals;
            else
                y(n+1, k+1) = (4^k * y(n+1, k) - y(n, k)) / (4^k - 1);
            end
        end

        % check if R(n, n) is close enough to R(n-1, n-1)
        if n > 0 && abs(y(n+1, n+1) - y(n, n)) < tol
            break;
        end
    end

    fprintf('--------------- romberg ---------------\n');
    fprintf('n: %d\n', n);
    fprintf('y = %0.10f\n', y(n+1));
    fprintf('function evaluations: %d\n', evals);
    fprintf('---------------------------------------\n');
end
