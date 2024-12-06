% 
% Jessie Li, CS 71 Fall 2023
%
% Exponential least squares approximation. 
% Fits y = b * e^(ax) to the data.
%
% Input:
%     x: x values
%     y: y values
%     linearized: if true, linearized, nonlinear otherwise
% 
% Returns:
%     coeffs: [b a]
%     err: sum of absolute errors

function [coeffs, err] = exponentialApprox(x, y, linearized)

    % linearized system: lny = lnb + ax
    lny = log(y);

    % coeffs = [lnb a]
    coeffs = polynomialApprox(x, lny, 1);

    % b = e^lnb
    coeffs(1) = exp(coeffs(1));

   if ~linearized
      % use linearized result as x0
      [~, coeffs, ~] = newtonNonlinear(@f, @J, coeffs, 20, 1e-10);
   end
    
   % calculate error
   yhat = coeffs(1) * exp(coeffs(2) * x);
   err = sum(abs(y - yhat));

    function z = f(t)
        c = t(2);
        d = t(1);

        z = [
            d * sum(exp(2*c*x)) - sum(y .* exp(c*x));
            d * sum(x .* exp(2*c*x)) - sum(x .* y .* exp(c*x))
         ];
    end
    
    function z = J(t)
       c = t(2);
       d = t(1);

       z = [    
           sum(exp(2*c*x)), 2 * d * sum(x .* exp(2*c*x)) - sum(x .* y .* exp(c*x));
           sum(x .* exp(2*c*x)), 2 * d * sum(x.^2 .* exp(2*c*x)) - sum(x.^2 .* y .* exp(c*x))
       ];
    end
end
