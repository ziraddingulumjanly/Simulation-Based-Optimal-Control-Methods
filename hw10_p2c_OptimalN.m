clear; clc; close all

% Problem constants
x0 = 1; tf = 1; rho = 2;

% Discretization sweep range
nt_range = 35:55;

% Reference analytical values (unconstrained)
J_star = -2/3;
u_star = 4/9;
x_star = 4/9;

fprintf('Constrained Optimal Control Problem\n');
fprintf('---------------------------------------------------------------\n');
fprintf('%5s | %10s | %9s | %9s | %9s\n', ...
    'n_t', 'J', 'u(t_f)', 'x(t_f)', 'Tot. Error');
fprintf('---------------------------------------------------------------\n');

for nt = nt_range
    auxdata.x0 = x0;
    auxdata.tf = tf;
    auxdata.rho = rho;
    auxdata.nt = nt;
    auxdata.t = linspace(0, auxdata.tf, auxdata.nt)';

    % Initial guess
    U0 = ones(auxdata.nt, 1) * 0.5;

    % Bounds
    lb = zeros(auxdata.nt,1);
    ub = ones(auxdata.nt,1) * 0.55;

    % Options
    options = optimoptions('fmincon', 'Display', 'off', 'MaxFunEvals', 1e4);

    % Optimize
    U_opt = fmincon(@(U) objective(U, auxdata), U0, [], [], [], [], ...
                    lb, ub, @(U) constraints(U, auxdata), options);

    % Simulate dynamics
    [t_sim, X_sim] = ode45(@(t,x) dynamics(t,x,U_opt,auxdata), ...
                           linspace(0, tf, 1000), x0);

    % Get final values
    final_J = objective(U_opt, auxdata);
    u_final = U_opt(end);
    x_final = X_sim(end);

    % Total error vs original analytical
    err = abs(final_J - J_star) + abs(u_final - u_star) + abs(x_final - x_star);

    % Print summary
    fprintf('%5d | %10.5f | %9.5f | %9.5f | %9.5f\n', ...
        nt, final_J, u_final, x_final, err);

    % Show full control and state 
    fprintf('   u(t): [');
    fprintf(' %.3f', U_opt(1:min(end,5)));  
    if nt > 5, fprintf(' ... %.3f', U_opt(end)); end
    fprintf(' ]\n');

    fprintf('   x(t_f): %.5f\n\n', x_final);
end

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
