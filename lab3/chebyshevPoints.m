% 
% Jessie Li, CS 71 Fall 2023
%
% Calculates the Chebyshev nodes on the interval [a, b]. 
%
% Parameters:
%     a: left bound of interval
%     b: right bound of interval
%     n: number of points

function p = chebyshevPoints(a, b, n)
    i = 0 : n-1; 
    xi = cos((2.*i + 1) .* pi / (2 .* (n+1))); 
    p = (a+b)/2 + (b-a)/2 .* xi;
end