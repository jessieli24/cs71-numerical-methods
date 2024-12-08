
%
% Jessie Li, CS 71 Fall 2023
% 
% Calculates the spectral radius of A.
%

function r = spectralRadius(A)
    n = size(A, 1);

    % T = -(D + L)^-1 * U
    T = -forwardSubstitute(tril(A), eye(n)) * triu(A, 1);
    r = max(abs(eig(T)));
end