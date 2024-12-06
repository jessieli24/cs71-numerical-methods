%
% Jessie Li, CS 71 Fall 2023
% 
% Approximates a definite integral with Gaussian quadrature.
% 
% Input:
%     f: function
%     a: lower limit
%     b: upper limit
%     n: number of points
%     ai: coefficients
%     xi: Gauss points
% 
% Returns:
%     y: approximated value
%     evals: number of function evaluations

function [y, evals] = gaussian(f, a, b, n, ai, xi)
    % map Gauss points from (-1, 1) to (a, b)
    yi = (a+b)/2 + (b-a)/2 .* xi;
        
    y = 0;
    evals = 0;

    for i = 1 : n
        y = y + ai(i) * f(yi(i));
        evals = evals + 1;
    end

    y = (b-a)/2 * y;
    
    fprintf('--------------- gaussian --------------\n');
    fprintf('y = %0.10f\n', y);
    fprintf('function evaluations: %d\n', evals);
    fprintf('---------------------------------------\n');
end