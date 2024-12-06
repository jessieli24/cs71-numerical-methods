% 
% Jessie Li, CS 71 Fall 2023
%
% Uses forward differences to approximate the second derivative. 
%
% Inputs:
%   y: y values
%   h: step size
%

function ddf = forwardDifference2(y, h)
    n = length(y) - 2;
    ddf = zeros(n, 1);

    for i = 1 : n
        ddf(i) = (y(i) - 2 * y(i+1) + y(i+2)) / h^2;
    end
end
