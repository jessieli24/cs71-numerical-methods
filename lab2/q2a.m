
% 
% q2a.m - visualizes function for Q2A
% 
% Jessie Li, CS 71 Fall 2023
%

% set default font to Times New Roman for all graphs
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

% --------------------- graph of f -------------------- %
% ----------------------------------------------------- %
x = linspace(0, 3, 500);
y = f(x);

% plot f
figure;
plot(x, y, 'LineWidth', 1.5);
xlabel('x');
ylabel('f(x)');
grid on;

saveas(gcf, 'results/q2a-f.png');
% -------------------- expected x --------------------- %
% ----------------------------------------------------- %
p = fzero(@f, [0 3]);
fprintf('x = %f\n', p);
fprintf('f''(x): %f\n', df(p));
fprintf('f''''(x) should be non-zero: %f\n', ddf(p));

function y = f(x)
y = (x + cos(x)) .* exp(-x.^2) + x.*cos(x);
end

function y = df(x)
y = cos(x) ...
    - exp(-x^2) * (sin(x) - 1) ...
    - x * sin(x) ...
    - 2 * x * exp(-x^2) * (x + cos(x));
end

function y = ddf(x)
y = 4 * x^2 * exp(-x^2) * (x + cos(x)) ...
    - 2 * exp(-x^2) * (x + cos(x)) ...
    - x * cos(x) ...
    - exp(-x^2) * cos(x) ...
    - 2 * sin(x) ...
    + 4 * x * exp(-x^2) * (sin(x) - 1);
end
