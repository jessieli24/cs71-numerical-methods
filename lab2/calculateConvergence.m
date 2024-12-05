% 
% Jessie Li, CS 71 Winter 2023
%
% Estimates the rate of convergence.
%
% a = [log(E[n+1])-log(E[n])] / [log(E[n])-log(E[n-1])]
%
% Parameters:
% - n: mumber of iterations (must be > 2)
% - err: error array
function a = calculateConvergence(n, err)
    if n < 3
        a = zeros(1, 1);
        return
    end 

    a = zeros(n-2, 1);
    logerr = log(err);

    prevdelta = logerr(2) - logerr(1);
    
    for i=2:n-1
       currdelta = logerr(i+1) - logerr(i);
       a(i-1) = currdelta / prevdelta;
       prevdelta = currdelta;
    end
end