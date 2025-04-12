clear; clc; close all;

% System parameters
alpha = 7.14;
beta  = 286.3;

A = [0    1;
     0 -alpha];
B = [0; beta];

% LQR weighting 
Q = [1 0; 
     0 0];

% Linear rho values
rho_vals = 1:10;

% Preallocate results
gain_list = zeros(length(rho_vals), 2);
pole_list = zeros(length(rho_vals), 2);

fprintf('LQR Results for rho = 1 to 10\n');
fprintf('-----------------------------------------------\n');
fprintf('%5s | %18s | %s\n', 'rho', 'K = [K1  K2]', 'Poles');
fprintf('-----------------------------------------------\n');

for i = 1:length(rho_vals)
    rho = rho_vals(i);
    R = rho;

    K = lqr(A, B, Q, R);
    Acl = A - B*K;
    poles = eig(Acl);

    % Store
    gain_list(i,:) = K;
    pole_list(i,:) = poles.';

    % Display
    fprintf('%5d | [%7.3f  %7.3f] | [%7.3f, %7.3f]\n', ...
        rho, K(1), K(2), real(poles(1)), real(poles(2)));
end

%  First Figure (Gains vs. rho) 
hf1 = figure('Position',[100 100 800 600]);  
hf1.Color = 'w'; 
hold on; grid off;

plot(rho_vals, gain_list(:,1), '-o', ...
    'LineWidth', 2, 'MarkerSize', 12, ...
    'MarkerFaceColor', [0 0.447 0.741], ...
    'DisplayName', '$K_1$');

plot(rho_vals, gain_list(:,2), '-o', ...
    'Color', [1 0 0], 'LineWidth', 2, ...
    'MarkerSize', 12, 'MarkerFaceColor', [1 0 0], ...
    'DisplayName', '$K_2$');

xlabel('$\rho$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('LQR Gain', 'Interpreter', 'latex', 'FontSize', 14);
title('LQR Gains $K_1$, $K_2$ vs. $\rho$', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');

set(gca, 'FontSize', 12, 'LineWidth', 1.2, ...
         'XColor', 'k', 'YColor', 'k');
box on
exportgraphics(hf1, 'p1a-gains-vs-rho.pdf', 'ContentType', 'vector');

% Keep system data in workspace
clearvars -except alpha beta A B Q

% Range of rho values on a log scale
rhoValues  = logspace(-3, 3, 100);
polesAll = zeros(2, length(rhoValues));

% Colors for each pole
c1 = [0.1 0.4 0.8];  
c2 = [0.2 0.7 0.2];  

% Loop over all rho values
for k = 1:length(rhoValues)
    Rk = rhoValues(k);
    Kk = lqr(A, B, Q, Rk);
    Acl_k = A - B * Kk;
    polesAll(:,k) = eig(Acl_k);
end

% Open-loop poles
polesOpenLoop = eig(A);

%  Second Figure (Full pole locus) 
hf2 = figure('Position',[150 150 800 600]);  
hf2.Color = 'w';
hold on;

% Plot open-loop poles
plot(real(polesOpenLoop), imag(polesOpenLoop), 'ks', ...
     'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Open-loop poles');

plot(real(polesAll(1,:)), imag(polesAll(1,:)), '-o', ...
     'Color', c1, 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', c1, 'DisplayName', 'Pole 1');
plot(real(polesAll(2,:)), imag(polesAll(2,:)), '-o', ...
     'Color', c2, 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', c2, 'DisplayName', 'Pole 2');

xlabel('Real Axis', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Imaginary Axis', 'Interpreter', 'latex', 'FontSize', 14);
title('Closed-Loop Pole Locus vs. $\rho$ (Full View)', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
axis equal
grid off
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'XColor', 'k', 'YColor', 'k');
box on
exportgraphics(hf2, 'p1b.pdf', 'ContentType', 'vector');

%  Third Figure (Zoomed-In View) 
hf3 = figure('Position',[200 200 800 600]);  
hf3.Color = 'w';
hold on;

% Plot open-loop poles
plot(real(polesOpenLoop), imag(polesOpenLoop), 'ks', ...
     'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Open-loop poles');

% Plot closed-loop poles
plot(real(polesAll(1,:)), imag(polesAll(1,:)), '-o', ...
     'Color', c1, 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', c1, 'DisplayName', 'Pole 1');
plot(real(polesAll(2,:)), imag(polesAll(2,:)), '-o', ...
     'Color', c2, 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', c2, 'DisplayName', 'Pole 2');

% Zoom in on x and y
xlim([-8 0.1])
ylim([-5 5])

xlabel('Real Axis', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Imaginary Axis', 'Interpreter', 'latex', 'FontSize', 14);
title('Closed-Loop Pole Locus vs. $\rho$ (Zoomed-In)', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
axis equal
grid off
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'XColor', 'k', 'YColor', 'k');
box on
exportgraphics(hf3, 'p1bzoomed.pdf', 'ContentType', 'vector');
