%
% Jessie Li, CS 71 Fall 2023
% 
% Solves Ax = b for a tridiagonal matrix A by performing Crout LU 
% decomposition on A. Assumes A is at least 2x2.
%

function x = solveTridiagonal(A, b)   
    A = LUDecomposeTridiagonal(A);

    % solve Ly = b using forward substitution
    x = forwardSubstitute(A, b);

    % solve Ux = y using backward substitution
    x = backwardSubstitute(A, x);
end
