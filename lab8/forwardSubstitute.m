%
% Jessie Li, CS 71 Fall 2023
% 
% Solves Lx = b with forward substitution.
%

function x = forwardSubstitute(A, b)
    for i = 1 : size(A, 1)
        b(i, :) = (b(i, :) - A(i, 1:i-1) * b(1:i-1, :)) / A(i, i);
    end

    x = b;
end