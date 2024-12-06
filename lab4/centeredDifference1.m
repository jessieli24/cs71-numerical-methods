% 
% Jessie Li, CS 71 Fall 2023
%
% Uses centered differences to approximate the first derivative. 
%
% Inputs:
%   y: y values
%   h: step size
%

function df = centeredDifference1(y, h)
    n = length(y) - 1;
    df = zeros(n-1, 1);
    
    for i = 2 : n
        df(i-1) = (y(i+1) - y(i-1)) / (2*h);
    end
end