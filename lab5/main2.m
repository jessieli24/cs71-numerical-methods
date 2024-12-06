% 
% main2.m - gaussian quadrature
% 
% Jessie Li, CS 71 Fall 2023
%

function main2(r_a01, r_a12, r_b01, r_b12)
    set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'DefaultTextInterpreter', 'latex');
    set(groot, 'DefaultLegendInterpreter', 'latex');

    % -------------------- constants ---------------------- %
    % ----------------------------------------------------- %
    MAX_STEPS = 11;
    
    % functions
    funcA = @(x) x.^2 .* exp(-x.^2);
    funcB = @(x) x.^(1/2) .* exp(-x.^2);

    a01 = zeros(MAX_STEPS, 1);
    a12 = zeros(MAX_STEPS, 1);
    b01 = zeros(MAX_STEPS, 1);
    b12 = zeros(MAX_STEPS, 1);

    % Gauss points and weights from:
    % https://pomax.github.io/bezierinfo/legendre-gauss.html
    
    % n = 1
    ai_xi1 = [2, 0];
    
    % n = 2
    ai_xi2 = [
        1.0000000000000000,	-0.5773502691896257;
        1.0000000000000000,	0.5773502691896257
    ];
    
    % n = 3
    ai_xi3 = [
        0.8888888888888888,	0.0000000000000000;
        0.5555555555555556,	-0.7745966692414834;
        0.5555555555555556,	0.7745966692414834
    ];
    
    % n = 4
    ai_xi4 = [
        0.6521451548625461,	-0.3399810435848563;
        0.6521451548625461,	0.3399810435848563;
        0.3478548451374538,	-0.8611363115940526;
        0.3478548451374538,	0.8611363115940526
    ];
    
    % n = 5
    ai_xi5 = [
        0.5688888888888889, 0.0000000000000000;
        0.4786286704993665,	-0.5384693101056831;
        0.4786286704993665,	0.5384693101056831;
        0.2369268850561891,	-0.9061798459386640;
        0.2369268850561891,	0.9061798459386640
    ];    
    
    % n = 6
    ai_xi6 = [
        0.3607615730481386,	0.6612093864662645;
        0.3607615730481386,	-0.6612093864662645;
	    0.4679139345726910,	-0.2386191860831969;
	    0.4679139345726910,	0.2386191860831969;
	    0.1713244923791704,	-0.9324695142031521;
	    0.1713244923791704,	0.9324695142031521
    ];
    
    % n = 7
    ai_xi7 = [
        0.4179591836734694,	0.0000000000000000;
	    0.3818300505051189,	0.4058451513773972;
	    0.3818300505051189,	-0.4058451513773972;
	    0.2797053914892766,	-0.7415311855993945;
	    0.2797053914892766,	0.7415311855993945;
	    0.1294849661688697,	-0.9491079123427585;
	    0.1294849661688697,	0.9491079123427585;
    ];
    
    % n = 8
    ai_xi8 = [
        0.3626837833783620,	-0.1834346424956498;
        0.3626837833783620,	0.1834346424956498;
        0.3137066458778873,	-0.5255324099163290;
        0.3137066458778873,	0.5255324099163290;
        0.2223810344533745,	-0.7966664774136267;
        0.2223810344533745,	0.7966664774136267;
        0.1012285362903763,	-0.9602898564975363;
        0.1012285362903763,	0.9602898564975363
    ];
    
    % n = 9
    ai_xi9 = [
        0.3302393550012598,	0.0000000000000000;
        0.1806481606948574,	-0.8360311073266358;
        0.1806481606948574,	0.8360311073266358;
        0.0812743883615744,	-0.9681602395076261;
        0.0812743883615744,	0.9681602395076261;
        0.3123470770400029,	-0.3242534234038089;
        0.3123470770400029,	0.3242534234038089;
        0.2606106964029354,	-0.6133714327005904;
        0.2606106964029354,	0.6133714327005904
    ];
    
    % n = 10
    ai_xi10 = [
        0.2955242247147529,	-0.1488743389816312;
        0.2955242247147529,	0.1488743389816312;
        0.2692667193099963,	-0.4333953941292472;
        0.2692667193099963,	0.4333953941292472;
        0.2190863625159820,	-0.6794095682990244;
        0.2190863625159820,	0.6794095682990244;
        0.1494513491505806,	-0.8650633666889845;
        0.1494513491505806,	0.8650633666889845;
        0.0666713443086881,	-0.9739065285171717;
        0.0666713443086881,	0.9739065285171717
    ];
    
    % n = 11
    ai_xi11 = [
        0.2729250867779006,	0.0000000000000000;
        0.2628045445102467,	-0.2695431559523450;
        0.2628045445102467,	0.2695431559523450;
        0.2331937645919905,	-0.5190961292068118;
        0.2331937645919905,	0.5190961292068118;
        0.1862902109277343,	-0.7301520055740494;
        0.1862902109277343,	0.7301520055740494;
        0.1255803694649046,	-0.8870625997680953;
        0.1255803694649046,	0.8870625997680953;
        0.0556685671161737,	-0.9782286581460570;
        0.0556685671161737,	0.9782286581460570;
    ];
    
    ai_xi = [
        ai_xi1; ai_xi2; ai_xi3; 
        ai_xi4; ai_xi5; ai_xi6; 
        ai_xi7; ai_xi8; ai_xi9; 
        ai_xi10; ai_xi11
    ];

    % ----------------------------------------------------- %
    % ----------------------------------------------------- %

    row = 1;
    for n = 1 : MAX_STEPS
        fprintf('---------------------- n = %d ----------------------\n', n);
        ai = ai_xi(row:row + n - 1, 1);
        xi = ai_xi(row:row + n - 1, 2);
     
        a01(n) = gaussian(funcA, 0, 1, n, ai, xi);
        a12(n) = gaussian(funcA, 1, 2, n, ai, xi);
        b01(n) = gaussian(funcB, 0, 1, n, ai, xi);
        b12(n) = gaussian(funcB, 1, 2, n, ai, xi);
        row = row + n;
    end
    
    table((1 : MAX_STEPS).', a01, a12, b01, b12)

    fprintf('function (a) with limits of integration from 0 to 1');
    plotGaussian(a01, MAX_STEPS, r_a01, 'results/a01')
    
    fprintf('function (a) with limits of integration from 1 to 2')
    plotGaussian(a12, MAX_STEPS, r_a12, 'results/a12')
    
    fprintf('function (b) with limits of integration from 0 to 1')
    plotGaussian(b01, MAX_STEPS, r_b01, 'results/b01')
    
    fprintf('function (b) with limits of integration from 1 to 2')
    plotGaussian(b12, MAX_STEPS, r_b12, 'results/b12')
end

% Plots results of the Gaussian quadrature.
% 
% Input:
%     yGaussian: Gaussian approximations
%     n: maximum number of points (1 to n)
%     yRomberg: Romberg approximations
%     filename: output file

function plotGaussian(yGaussian, n, yRomberg, filename)
    figure
    plot(1:n, yGaussian, '-o')
    yline(yRomberg, '--')   % yline for Romberg

    xlabel('n')
    ylabel('Approximated integral value')
    title('Gaussian Integral Approximation v. Number of Points')

    % saveas(gcf, sprintf('%s-g-approx.png', filename))

    % plot differences v. n
    delta = zeros(n-1);

    for i = 1:n-1
        delta(i) = abs(yGaussian(i+1) - yGaussian(i));
    end

    figure
    plot(2:n, delta, '-o')
    yline(1e-9, '--')

    set(gca, 'YScale', 'log')

    xlabel('n')
    ylabel('Absolute difference')
    title({'Absolute Difference Between Successive', 'Integral Approximations v. n'})
    
    % saveas(gcf, sprintf('%s-g-diff.png', filename))
end
