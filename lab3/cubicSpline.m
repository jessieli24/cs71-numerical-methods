% 
% Jessie Li, CS 71 Fall 2023
%
% Natural cubic spline. 
% Based on Algorithm 3.4 in the textbook.
%
% Parameters:
%     n: number of points
%     x: x values as a list of size n
%     y: y values as a list of size n
% 
% Returns:
%     [a, b, c, d]: four (n-1) x 1 arrays of polynomial coefficients:
%                   a + b(x - xi) + c(x - xi)^2 + d(x - xi)^3

function [a, b, c, d] = cubicSpline(n, x, y)
    a = y;
    m = n - 1;    % number of segments
    
    h = zeros(m, 1);
    for i = 1 : m
        h(i) = x(i+1) - x(i);
    end
    
    p = zeros(m, 1);
    for i = 2 : m
        p(i) = 3/h(i) * (a(i+1) - a(i)) - 3/h(i-1) * (a(i) - a(i-1));
    end
    
    l = zeros(n, 1);
    u = zeros(n, 1);
    z = zeros(n, 1);
    
    l(1) = 1;
    l(n) = 1;
    
    for i = 2 : m
        l(i) = 2 * (x(i+1) - x(i-1)) - h(i-1) * u(i-1);
        u(i) = h(i) / l(i);
        z(i) = (p(i) - h(i-1) * z(i-1)) / l(i);
    end
    
    b = zeros(n, 1);
    c = zeros(n, 1);
    d = zeros(n, 1);
    
    for j = m : -1 : 1
        c(j) = z(j) - u(j)*c(j+1);
        b(j) = (a(j+1) - a(j)) / h(j) - h(j) * (c(j+1) + 2*c(j)) / 3;
        d(j) = (c(j+1) - c(j)) / (3 * h(j));
    end
    
     a = a(1:end-1);
     b = b(1:end-1);
     c = c(1:end-1);
     d = d(1:end-1);
end