clear; clc; close all;

% Parameters
x0 = 1; tf = 1; rho = 2;
nt_range = 3:35;

% True optimal values
J_star = -2/3;
u_star = 4/9;
x_star = 4/9;

results = [];

fprintf('Optimal Control Simulation (Sweep over n_t)\n');
fprintf('-------------------------------------------------------------\n');
fprintf('%5s | %10s | %9s | %9s | %9s\n', ...
    'n_t', 'J', 'u(t_f)', 'x(t_f)', 'Tot. Error');
fprintf('-------------------------------------------------------------\n');

for nt = nt_range
    auxdata.x0 = x0;
    auxdata.tf = tf;
    auxdata.rho = rho;
    auxdata.nt = nt;
    auxdata.t = linspace(0, auxdata.tf, auxdata.nt)';

    U0 = ones(auxdata.nt, 1) * 0.5;
    options = optimoptions('fmincon', 'Display', 'off', 'MaxFunEvals', 1e4);

    U_opt = fmincon(@(U) objective(U, auxdata), U0, [], [], [], [], ...
                    zeros(auxdata.nt,1), ones(auxdata.nt,1), ...
                    @(U) constraints(U, auxdata), options);

    final_J = objective(U_opt, auxdata);
    [~, X_sim] = ode45(@(t,x) dynamics(t,x,U_opt,auxdata), ...
                       linspace(0, tf, 1000), x0);

    u_final = U_opt(end);
    x_final = X_sim(end);

    % Errors
    errJ = abs(final_J - J_star);
    errU = abs(u_final - u_star);
    errX = abs(x_final - x_star);
    total_err = errJ + errU + errX;

    results = [results; nt, final_J, u_final, x_final, total_err];

    % Print
    fprintf('%5d | %10.5f | %9.5f | %9.5f | %9.5f\n', ...
        nt, final_J, u_final, x_final, total_err);
end

% Extract
nt_vals = results(:,1);
J_vals = results(:,2);
u_vals = results(:,3);
x_vals = results(:,4);
err_vals = results(:,5);

% Colors
c1 = [0 0.447 0.741]; % Blue
c2 = [0.850 0.325 0.098]; % Orange
c3 = [0.466 0.674 0.188]; % Green
c4 = [0.494 0.184 0.556]; % Purple

% Plotting
hf = figure('Position', [100 100 1200 700]); 
hf.Color = 'w';
tiledlayout(2,2)

%  J 
nexttile
plot(nt_vals, J_vals, '-o', 'Color', c1, ...
     'MarkerFaceColor', c1, 'LineWidth', 2)
yline(J_star, '--r', '$J^*$', 'Interpreter', 'latex');
title('Objective $J$ vs $n_t$', 'Interpreter', 'latex')
xlabel('$n_t$', 'Interpreter', 'latex'); ylabel('$J$', 'Interpreter', 'latex')
grid on

%  u(t_f) 
nexttile
plot(nt_vals, u_vals, '-o', 'Color', c2, ...
     'MarkerFaceColor', c2, 'LineWidth', 2)
yline(u_star, '--r', '$u^*(1)$', 'Interpreter', 'latex');
title('Final Control $u(t_f)$ vs $n_t$', 'Interpreter', 'latex')
xlabel('$n_t$', 'Interpreter', 'latex'); ylabel('$u(t_f)$', 'Interpreter', 'latex')
grid on

%  x(t_f) 
nexttile
plot(nt_vals, x_vals, '-o', 'Color', c3, ...
     'MarkerFaceColor', c3, 'LineWidth', 2)
yline(x_star, '--r', '$x^*(1)$', 'Interpreter', 'latex');
title('Final State $x(t_f)$ vs $n_t$', 'Interpreter', 'latex')
xlabel('$n_t$', 'Interpreter', 'latex'); ylabel('$x(t_f)$', 'Interpreter', 'latex')
grid on

%  Error 
nexttile
plot(nt_vals, err_vals, '-s', 'Color', c4, ...
     'MarkerFaceColor', c4, 'LineWidth', 2)
title('Combined Error vs $n_t$', 'Interpreter', 'latex')
xlabel('$n_t$', 'Interpreter', 'latex'); ylabel('Total Absolute Error', 'Interpreter', 'latex')
grid on
exportgraphics(hf, 'mine_problem_convergence.pdf', 'ContentType', 'vector');

% Functions
function J = objective(U, auxdata)
    [~, X] = ode45(@(t,x) dynamics(t,x,U,auxdata), auxdata.t, auxdata.x0);
    x_vals = X(:,1);
    x_vals(x_vals < 1e-6) = 1e-6;
    L = (U.^2)./x_vals - auxdata.rho * U;
    J = trapz(auxdata.t, L);
end

function dx = dynamics(t, x, U, auxdata)
    u = interp1(auxdata.t, U, t, 'linear', 'extrap');
    dx = -u;
end

function [g, h] = constraints(U, auxdata)
    [~, X] = ode45(@(t,x) dynamics(t,x,U,auxdata), auxdata.t, auxdata.x0);
    x_vals = X(:,1);
    epsilon = 1e-6;
    g = [-x_vals + epsilon; x_vals - auxdata.x0];
    h = [];
end
