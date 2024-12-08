%
% Jessie Li, CS 71 Fall 2023
% 
% Solves Ax = b for a tridiagonal matrix A by performing Crout LU 
% decomposition on A. Assumes A is at least 2x2.
%

function x = solveTridiagonal(A, b)   
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

    % solve Ly = b using forward substitution
    for i = 1 : n
        for j = 1 : i-1
            b(i) = b(i) - A(i, j) * b(j);
        end
        
        b(i) = b(i) / A(i, i);
    end

    % solve Ux = y using backward substitution
    for i = n : -1 : 1
        for j = i+1 : n
            b(i) = b(i) - A(i, j) * b(j);
        end
        
    end

    x = b;
end
