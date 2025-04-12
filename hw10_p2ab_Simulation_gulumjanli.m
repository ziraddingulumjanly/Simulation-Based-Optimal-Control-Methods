clear; clc; close all

% Problem parameters
auxdata.x0 = 1;     
auxdata.tf = 1;      
auxdata.rho = 2;     
auxdata.nt = 11;     % Time discretization points
auxdata.t = linspace(0, auxdata.tf, auxdata.nt)';

U0 = ones(auxdata.nt, 1) * 0.5;
options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunEvals', 1e4);

% Solve optimal control problem
fprintf('Solving optimal control problem...\n')
U_opt = fmincon(@(U) objective(U, auxdata), U0, [], [], [], [], ...
                zeros(auxdata.nt,1), ones(auxdata.nt,1), ...
                @(U) constraints(U, auxdata), options);

% Simulate optimal trajectory
[t_sim, X_sim] = ode45(@(t,x) dynamics(t,x,U_opt,auxdata), ...
                       linspace(0, auxdata.tf, 1000), auxdata.x0);

% Plot state x(t)
c1 = [0.1 0.4 0.8];  

hf1 = figure; hold on; hf1.Color = 'w';
plot(t_sim, X_sim, '-', 'Color', c1, 'LineWidth', 2)
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$x(t)$ (state)', 'Interpreter', 'latex', 'FontSize', 14)
%title('State Evolution Over Time', 'Interpreter', 'latex', 'FontSize', 16)
grid off
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'XColor', 'k', 'YColor', 'k')
exportgraphics(hf1, 'p2state.pdf', 'ContentType', 'vector')

% Plot optimal control u(t)
c2 = [0.8 0.1 0.1];  

hf2 = figure; hold on; hf2.Color = 'w';
plot(auxdata.t, U_opt, 'o-', 'Color', c2, ...
    'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', c2)
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$u(t)$ (control)', 'Interpreter', 'latex', 'FontSize', 14)
%title('Optimal Control Input Over Time', 'Interpreter', 'latex', 'FontSize', 16)
grid off
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'XColor', 'k', 'YColor', 'k')
exportgraphics(hf2, 'p2control.pdf', 'ContentType', 'vector')

% Display results
total_control = trapz(auxdata.t, U_opt);
fprintf('\nTotal integrated control: %.4f\n', total_control)

% Objective Function
function J = objective(U, auxdata)
    [~, X] = ode45(@(t,x) dynamics(t,x,U,auxdata), auxdata.t, auxdata.x0);
    x_vals = X(:,1);
    x_vals(x_vals < 1e-6) = 1e-6;  
    L = (U.^2)./x_vals - auxdata.rho * U;
    J = trapz(auxdata.t, L);
end

% Dynamics Function
function dx = dynamics(t, x, U, auxdata)
    u = interp1(auxdata.t, U, t, 'linear', 'extrap');
    dx = -u;
end

% Path Constraints
function [g, h] = constraints(U, auxdata)
    [~, X] = ode45(@(t,x) dynamics(t,x,U,auxdata), auxdata.t, auxdata.x0);
    x_vals = X(:,1);
    epsilon = 1e-6;
    g = [-x_vals + epsilon; x_vals - auxdata.x0]; 
    h = [];
end
