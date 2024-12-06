% 
% Jessie Li, CS 71 Fall 2023
%
% Polynomial least squares approximation. 
%
% Input:
%   x: x data values
%   y: y data values
%   n: order of approximating polynomial
%
% Returns:
%   coeffs: coefficients of the least squares polynomial
%       p_n(x) = a_0 + a_1 * x + ... + a_n * x^n
%   err: sum of absolute errors

function [coeffs, err] = polynomialApprox(x, y, n)
    X = zeros(n+1, n+1);
    b = zeros(n+1, 1);

    for k = 0 : n
        for j = 0 : n
            X(k + 1, j + 1) = sum(x.^(k + j));        
        end
        
        b(k + 1) = sum(y .* x.^k);
    end  

    % solve for coefficients
    coeffs = X \ b;

    % calculate error
    yhat = polyval(flip(coeffs, 1), x);
    err = sum(abs(y - yhat));
end