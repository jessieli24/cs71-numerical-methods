% 
% main1.m - polynomial and exponential approximations
% 
% Jessie Li, CS 71 Fall 2023
%

set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

% data points
data = [
    0, 2.5237305e+03;
    5, 2.0864229e+02;
    6, 1.8697888e+02;
    13, 3.9214807e+00;
    16, 1.5271976e+00;
    22, 6.8035297e-02;
    35, 7.3126314e-05;
    38, 1.4585192e-05;
    42, 3.0299061e-06;
    44, 8.6349502e-07;
    47, 2.4089013e-07;
    50, 4.6306127e-08
];

% separate x and y values
x = data(:, 1);
y = data(:, 2);

% linear polynomial fit
[coeffs, error] = polynomialApprox(x, y, 1)
plotPolynomial(x, y, coeffs, 'Linear Polynomial Fit');

% cubic polynomial fit
[coeffs, error] = polynomialApprox(x, y, 3)
plotPolynomial(x, y, coeffs, 'Cubic Polynomial Fit');

% linear polynomial fit to linearized data
[coeffs, error] = exponentialApprox(x, y, true)
plotExponential(x, y, coeffs, 'Linear Fit to Linearized Data');

% nonlinear exponential fit
[coeffs, error] = exponentialApprox(x, y, false)
plotExponential(x, y, coeffs, 'Nonlinear Fit');

function f = plotPolynomial(x, y, coeffs, graphTitle)
    xx = min(x) : 1e-3 : max(x);
    yy = polyval(flip(coeffs, 1), xx);
    
    f = figure;
    plot(x, y, 'o', xx, yy);
    
    xlabel('x');
    ylabel('y');
    title(graphTitle);
end

function f = plotExponential(x, y, coeffs, graphTitle)
    a = coeffs(2);
    b = coeffs(1);

    xx = min(x) : 1e-3 : max(x);
    yy = b .* exp(a * xx);
    
    f = figure;
    plot(x, y, 'o', xx, yy);
    set(gca, 'YScale', 'log');
    
    xlabel('x');
    ylabel('y');
    title(graphTitle);
end
