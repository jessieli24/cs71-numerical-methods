%
% Jessie Li, CS 71 Fall 2023
% 
% Solves Ax = b using Gauss-Seidel iteration.
%
% Inputs:
%   A: matrix
%   b: vector
%   T0: initial value
%   tol: tolerance
%   maxIters: maximum number of iterations
% 
% Returns:
%   xCurr: approximate solution x
%   r: spectral radius estimate

function [xCurr, r] = gaussSeidel(A, b, T0, tol, maxIters)
    xPrev = T0;
    xCurr = T0;
    deltaPrev = 0;
    r = 0;

    for k = 1 : maxIters
        for i = 1 : size(xCurr, 1)
            xCurr(i) = (b(i) - A(i, 1:i-1) * xCurr(1:i-1) - A(i, i+1:end) * xCurr(i+1:end)) / A(i, i);
        end

        % L∞ relative norm stopping criterion 
        error = max(abs((xCurr - xPrev) ./ xCurr));
        deltaCurr = abs(xCurr - xPrev);

        if error < tol
            fprintf('n + 1 = %d: Gauss-Seidel converged after %d iterations.\n', size(A, 1)+1, k);
            % estimate spectral radius 
            r = estimateSpectralRadius(deltaPrev, deltaCurr);
            break;

        elseif k == maxIters
            fprintf('n + 1 = %d: Gauss-Seidel reached maximum of %d iterations, error = %.5f\n', size(A, 1)+1, maxIters, error);
            r = estimate_spectral_radius(deltaPrev, deltaCurr);

        end

        xPrev = xCurr;
        deltaPrev = deltaCurr;
    end
end
