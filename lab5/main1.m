% 
% main1.m - numerical integration
% 
% Jessie Li, CS 71 Fall 2023
%

function [r_a01, r_a12, r_b01, r_b12] = main1()
    set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'DefaultTextInterpreter', 'latex');
    set(groot, 'DefaultLegendInterpreter', 'latex');
    
    % -------------------- constants ---------------------- %
    % ----------------------------------------------------- %
    TOL = 1e-9;
    MAX_STEPS = 25;
    
    % functions
    funcA = @(x) x.^2 .* exp(-x.^2);
    funcB = @(x) x.^(1/2) .* exp(-x.^2);
    % ----------------------------------------------------- %
    % ----------------------------------------------------- %
    
    % plot functions
    x = 0:1e-3:2;
    
    figure
    colororder('reef')
    plot(x, funcA(x), x, funcB(x))
    xline(1, '--')
    
    legend('a', 'b', 'Location', 'northeast')
    
    xlabel('x')
    ylabel('y')
    
    % saveas(gcf, 'results/q1-functions.png')
    
    % --------------------- function a -------------------- %
    % ----------------------------------------------------- %
    
    % function (a) with limits of integration from 0 to 1
    fprintf('--------------------------------------------------\n')
    fprintf('function (a) with limits of integration from 0 to 1\n')
    
    [yr, nr] = romberg(funcA, 0, 1, TOL, MAX_STEPS);
    [yt, nt] = evaluateTrapezoidal(funcA, 0, 1, TOL, MAX_STEPS);
    
    plotRTIntegration(yr, nr, yt, nt, 'results/q1-a01')
    
    r_a01 = yr(nr+1, nr+1);
    
    % function (a) with limits of integration from 1 to 2
    fprintf('--------------------------------------------------\n')
    fprintf('function (a) with limits of integration from 1 to 2\n')
    
    [yr, nr] = romberg(funcA, 1, 2, TOL, MAX_STEPS);
    [yt, nt] = evaluateTrapezoidal(funcA, 1, 2, TOL, MAX_STEPS);
    
    plotRTIntegration(yr, nr, yt, nt, 'results/q1-a12')
    
    r_a12 = yr(nr+1, nr+1);
    
    % --------------------- function b -------------------- %
    % ----------------------------------------------------- %
    
    % function (b) with limits of integration from 0 to 1
    fprintf('--------------------------------------------------\n')
    fprintf('function (b) with limits of integration from 0 to 1\n')
    
    [yr, nr] = romberg(funcB, 0, 1, TOL, MAX_STEPS);
    [yt, nt] = evaluateTrapezoidal(funcB, 0, 1, TOL, MAX_STEPS);
    
    plotRTIntegration(yr, nr, yt, nt, 'results/q1-b01')
    
    r_b01 = yr(nr+1, nr+1);
    
    % function (b) with limits of integration from 1 to 2
    fprintf('--------------------------------------------------\n')
    fprintf('function (b) with limits of integration from 1 to 2\n')
    
    [yr, nr] = romberg(funcB, 1, 2, TOL, MAX_STEPS);
    [yt, nt] = evaluateTrapezoidal(funcB, 1, 2, TOL, MAX_STEPS);
    
    plotRTIntegration(yr, nr, yt, nt, 'results/q1-b12')
    
    r_b12 = yr(nr+1, nr+1);
end

% ---------------- helper functions ------------------- %
% ----------------------------------------------------- %

% Plots integral approximations v. n. and absolute differences v. n.
% 
% Input:
%     rtable: Romberg table
%     nr: number of steps (2^n)
%     ttable: trapezoidal approximations
%     nt: number of steps (2^n)
%     filename: output file

function plotRTIntegration(rtable, nr, ttable, nt, filename)
    % get Romberg approximations along the diagonal
    yr = zeros(nr+1, 1);
    for i = 1 : nr+1
        yr(i) = rtable(i, i);
    end

    yt = ttable(1:nt + 1);

    % plot approximations v. n
    figure

    hold on
    plot(0:nr, yr, '-o')
    plot(0:nt, yt, '-o')
    hold off
    
    legend('romberg', 'trapezoidal', 'Location', 'northeast')

    xlabel('n')
    ylabel('Approximate integral value')
    title({'Approximate Integral Value', 'v. Number of Steps'})
    
    % saveas(gcf, sprintf('%s-rt-approx.png', file))

    % compute errors
    rErr = zeros(nr, 1);
    for i = 1 : nr
        rErr(i) = abs(rtable(i+1, i+1) - rtable(i, i));
    end

    tErr = zeros(nt, 1);
    for i = 1 : nt
        tErr(i) = abs(ttable(i+1) - ttable(i));
    end

    % plot errors v. n
    figure

    hold on
    plot(1:nr, rErr, '-o')
    plot(1:nt, tErr, '-o')
    yline(1e-9, '--')
    hold off
    
    set(gca, 'YScale', 'log')
    legend('romberg', 'trapezoidal', 'Location', 'northeast')

    xlabel('n')
    ylabel('Absolute difference')
    title({'Absolute Difference Between Successive', 'Integral Approximations v. n'})

    % saveas(gcf, sprintf('%s-rt-diff.png', file))
end