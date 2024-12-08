%
% Jessie Li, CS 71 Fall 2023
% 
% Calculates the matrix for solving the bioheat BVP with the 
% Crank-Nicholson time-stepping algorithm. 
%

function [A, B, c, h] = getCrankNicholsonMatrix(n, dt)
    L = 1;
    lambda2 = 2.7;
    Ta = 37;
    Tc = 37;
    Ts = 32;

    A = zeros(n, n);
    B = zeros(n, n);

    c = zeros(n, 1);
    h = L / (n + 1);

    % Crank-Nicholson time-stepping algorithm
    for i = 1 : n
        A(i, i) = 2*h^2 + 2*dt + dt*h^2*lambda2;
        B(i, i) = 2*h^2 - 2*dt - dt*h^2*lambda2;
    
        if i < n
            A(i, i+1) = -dt;
            B(i, i+1) = dt;
        end

        if i > 1
            A(i, i-1) = -dt;
            B(i, i-1) = dt;
        end
    end

    % boundary conditions
    c(1) = (Tc - Ta) * 2 * dt;
    c(n) = (Ts - Ta) * 2 * dt;
end
