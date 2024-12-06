% 
% Jessie Li, CS 71 Fall 2023
%
% Uses centered differences to approximate the second derivative. 
%
% Inputs:
%   y: y values
%   h: step size
%

function ddf = centeredDifference2(y, h)
    n = length(y) - 1;
    ddf = zeros(n-1, 1);
    
    for j = 2 : n
        ddf(j-1) = (y(j-1) - 2 * y(j) + y(j+1)) / h^2;
    end
end