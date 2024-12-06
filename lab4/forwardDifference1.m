% 
% Jessie Li, CS 71 Fall 2023
%
% Uses forward differences to approximate the first derivative. 
%
% Inputs:
%   y: y values
%   h: step size
%

function df = forwardDifference1(y, h)
    n = length(y) - 1;
    df = zeros(n, 1);
    
    for i = 1 : n
        df(i) = (y(i+1) - y(i)) / h;
    end
end
