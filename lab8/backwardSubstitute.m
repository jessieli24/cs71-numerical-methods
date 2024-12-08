%
% Jessie Li, CS 71 Fall 2023
% 
% Solves Ux = b with backward substitution.
% Assumes diagonal entries A(i, i) = 1.
%

function x = backwardSubstitute(A, b)
    n = size(A, 1);
    
    for i = n : -1 : 1
        b(i, :) = (b(i, :) - A(i, i+1:n) * b(i+1:n, :));
    end

    x = b;
end