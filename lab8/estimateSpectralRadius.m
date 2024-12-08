%
% Jessie Li, CS 71 Fall 2023
% 
% Estimates the spectral radius of a matrix A.
%

function r = estimateSpectralRadius(deltaPrev, deltaCurr)
    deltaPrevNorm = sqrt(sum(abs(deltaPrev .^ 2)));
    deltaCurrNorm = sqrt(sum(abs(deltaCurr .^ 2)));

    if deltaPrevNorm ~= 0
        r = deltaCurrNorm / deltaPrevNorm;
    else
        r = 0;
    end
end