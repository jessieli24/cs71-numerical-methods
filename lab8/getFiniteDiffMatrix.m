%
% Jessie Li, CS 71 Fall 2023
% 
% Gets the finite difference matrix representing a system of finite 
% difference equations Ax = b with n interior nodes, for the bioheat BVP 
% describing the steady-state temperature distribution in the absence of 
% heating.
%

function [A, b, h] = getFiniteDiffMatrix(n)
    L = 1;
    lambda2 = 2.7;
    Ta = 37;
    Tc = 37;
    Ts = 32;

    A = zeros(n, n);
    b = zeros(n, 1);
    h = L / (n + 1);

    for i = 1 : n
        A(i, i) = -(2 + h^2 * lambda2);

        if i < n
            A(i, i+1) = 1;
        end

        if i > 1
            A(i, i-1) = 1;
        end
    end

    % boundary conditions
    b(1) = -(Tc - Ta);
    b(n) = -(Ts - Ta);
end