%
% Jessie Li, CS 71 Fall 2023
% 
% Wrapper for compositeTrapezoidal.m.
% 
% Input:
%     f: function
%     a: lower limit
%     b: upper limit
%     tol: tolerance
%     max_steps: maximum number of iterations
% 
% Returns:
%     y: array of approximations, indexed by step number
%     n: number of times subintervals were halved
%     evals: number of function evaluations

function [y, n, evals] = evaluateTrapezoidal(f, a, b, tol, max_steps)
    y = zeros(max_steps+1, 1);
    evals = 0;
    
    for n = 0 : max_steps
        [y(n+1), evals] = compositeTrapezoidal(f, a, b, 2^n);
        
        % Check if y(n) is close enough to y(n-1)
        if n > 0 && abs(y(n+1) - y(n)) < tol
            break;
        end
    end

    fprintf('------------- trapezoidal -------------\n');
    fprintf('n: %d\n', n);
    fprintf('y = %0.10f\n', y(n+1));
    fprintf('function evaluations: %d\n', evals);
    fprintf('---------------------------------------\n');
end