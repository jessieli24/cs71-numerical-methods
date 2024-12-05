% 
% Jessie Li, CS 71 Fall 2023
%
% Neville's method. 
%
% Parameters:
%   xin: input x data values as an array of size n
%   yin: input y data values as an array of size n
%   n: number of points (order n-1)
%   x: x value to interpolate

function y = neville(xin, yin, n, x)
    
    p = yin; % first row

    for i = 1 : n-1
        for j = 1 : n-i
            p(j) = ((x - xin(j+i)) * p(j) - ...
                    (x - xin(j)) * p(j+1)) / ...
                    (xin(j) - xin(j+i));
        end
    end
    
    y = p(1);
end