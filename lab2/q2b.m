
% 
% q2b.m - visualizes function for Q2B 
% 
% Jessie Li, CS 71 Fall 2023
%

% set default font to Times New Roman for all graphs
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

% --------------------- graph of f -------------------- %
% ----------------------------------------------------- %
x = linspace(0, 3, 200);
y = f(x);

% plot f
figure;
plot(x, y, 'LineWidth', 1.5);
xlabel('x');
ylabel('f(x)');
grid on;

saveas(gcf, 'results/q2b-f.png');

% -------------------- expected x --------------------- %
% ----------------------------------------------------- %
% p = fzero(@f, [1.5 2]); % should be the same as 2A
p = 1.636723;
fprintf('x = %f (from 2A)\n', p);
fprintf('f''(x) should be zero: %f\n', df(p));
fprintf('f''''(x) should be non-zero: %f\n', ddf(p));

function y = f(x)
y = ((x + cos(x)) .* exp(-x.^2) + x.*cos(x)).^2;
end

function y = df(x)
y = -2 * (exp(-x^2) * (x + cos(x)) + x * cos(x)) * ...
    (exp(-x^2) * (sin(x) - 1) - cos(x) + x * sin(x) + ...
    2 * x * exp(-x^2) * (x + cos(x)));
end

function y = ddf(x)
y = 2 * (exp(-x^2) * (sin(x) - 1) - cos(x) + x * sin(x) + ...
    2 * x * exp(-x^2) * (x + cos(x)))^2 - 2 * (exp(-x^2) * ...
    (x + cos(x)) + x * cos(x)) * (2 * sin(x) + 2 * exp(-x^2) * ...
    (x + cos(x)) + x * cos(x) + exp(-x^2) * cos(x) - 4 * x^2 * ...
    exp(-x^2) * (x + cos(x)) - 4 * x * exp(-x^2) * (sin(x) - 1));
end
