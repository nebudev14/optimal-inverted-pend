clc;
clear;
addpath('~/packages/casadi-matlab/')
import casadi.*

N = 20; % number of interior points
opti = casadi.Opti();

g = 10; % gravity
l = 4; % meters
m = 2; % kg

X = opti.variable(2, N+1); % State space trajectory
theta = X(1, :);
theta_dot = X(2, :);
U = opti.variable(1, N); % control inputs
T = 5.0; % time
J = 0;
for k=1:N
    J = J + U(k)^2;
end

opti.minimize(J);

f = @(x, u) [x(2); ((g/l)*sin(x(1))) + (u/(m*l^2))]; % dynamics of system
dt = T/N;
% multiple shooting via RK4
for k=1:N
   k1 = f(X(:,k),         U(:,k));
   k2 = f(X(:,k)+dt/2*k1, U(:,k));
   k3 = f(X(:,k)+dt/2*k2, U(:,k));
   k4 = f(X(:,k)+dt*k3,   U(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

% set some arbitrary control constraint
% opti.subject_to(0 <= U <= 2);

opti.subject_to(theta(1) == 0);
opti.subject_to(theta_dot(1) == 0);
opti.subject_to(theta(N+1) == pi);
opti.subject_to(U(1) == 0);
% 
% opti.subject_to(T>=0);

opti.set_initial(theta_dot, 2);
opti.set_initial(U, 1);

opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solve

t_grid = linspace(0, T, N+1);

figure;
subplot(3,1,1);
plot(t_grid, sol.value(theta), '-o');
xlabel('Time [s]'); ylabel('\theta [rad]');
title('Optimal Angle Trajectory');

subplot(3,1,2);
plot(t_grid, sol.value(theta_dot), '-o');
xlabel('Time [s]'); ylabel('d\theta/dt [rad/s]');
title('Optimal Angular Velocity Trajectory');

subplot(3,1,3);
plot(linspace(0,T,N), sol.value(U), '-o');
xlabel('Time [s]'); ylabel('u [N*m]');
title('Optimal Control Input');