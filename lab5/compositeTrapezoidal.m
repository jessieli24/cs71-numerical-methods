%
% Jessie Li, CS 71 Fall 2023
% 
% Approximates a definite integral with the composite trapezoidal rule.
% 
% Input:
%     f: function
%     a: lower limit
%     b: upper limit
%     n: number of subintervals
% 
% Returns:
%     y:  approximated value of integral
%     evals: number of function evaluations

function [y, evals] = compositeTrapezoidal(f, a, b, n)
    h = (b - a)/n;
    evals = n + 1;
    y = a : h : b; 
    
    y(2:n) = 2 * f(y(2:end-1)); % interior points
    y(1) = f(y(1));             % left endpoint
    y(n+1) = f(y(n+1));         % right endpoint

    y = sum(y) * h/2;
end
