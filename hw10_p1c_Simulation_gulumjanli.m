clear; clc; close all

%  System parameters 
alpha = 7.14;
beta = 286.3;

A = [0 1;
     0 -alpha];
B = [0;
     beta];
C = [1 0];
Q = [1 0;
     0 0];

rho = 0.0231013;      % Selected rho
R = rho;
x0 = [0.3; 0];        

%  LQR gain and closed-loop system 
K = lqr(A, B, Q, R);
Acl = A - B*K;

%  Simulation time 
t = linspace(0, 0.2, 1000);  

%  Simulate state 
[~, X] = ode45(@(t, x) Acl * x, t, x0);
theta = X(:,1);           
u = -K * X.';             
u = u.';                  

%  Plot theta(t) 
hf1 = figure; hf1.Color = 'w'; hold on
plot(t * 1000, theta, 'r-', 'LineWidth', 2)
% xline(50, '--k', '50 ms', 'LabelVerticalAlignment', 'bottom', ...
%        'Interpreter','latex', 'FontSize', 11)
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\theta(t)$ [rad]', 'Interpreter', 'latex', 'FontSize', 14)
%title('$\theta(t)$ Response with LQR ($\rho = 0.5$)', 'Interpreter', 'latex')
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'XColor', 'k', 'YColor', 'k')
grid off; box on
exportgraphics(hf1, 'p1c_theta.pdf', 'ContentType', 'vector');

%  Plot u(t) 
hf2 = figure; hf2.Color = 'w'; hold on
plot(t * 1000, u, 'b-', 'LineWidth', 2)
% xline(50, '--k', '50 ms', 'LabelVerticalAlignment', 'bottom', ...
%        'Interpreter','latex', 'FontSize', 11)
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$u(t)$ [V]', 'Interpreter', 'latex', 'FontSize', 14)
%title('Control Input $u(t)$ for LQR ($\rho = 0.5$)', 'Interpreter', 'latex')
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'XColor', 'k', 'YColor', 'k')
grid off; box on
exportgraphics(hf2, 'p1c_control.pdf', 'ContentType', 'vector');
