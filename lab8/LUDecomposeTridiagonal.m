%
% Jessie Li, CS 71 Fall 2023
% 
% Performs Crout LU decomposition on a tridiagonal matrix A.
%

function A = LUDecomposeTridiagonal(A)
    n = size(A, 1);

    % U(1, 2) = A(1, 2) / L(1, 1)
    A(1, 2) = A(1, 2)/A(1, 1);

    for i = 2 : n-1
        % A(i, i-1) = L(i, i-1)

        % L(i, i) = A(i, i) - L(i, i-1) * U(i-1, i)
        A(i, i) = A(i, i) - A(i, i-1) * A(i-1, i); 

        % U(i, i+1) = A(i, i+1)/ L(i, i)
        A(i, i+1) = A(i, i+1) / A(i, i);
    end

    % L(n, n) = A(n, n) - L(n, n-1) * U(n-1, n)
    A(n, n) = A(n, n) - A(n, n-1) * A(n-1, n);
end
