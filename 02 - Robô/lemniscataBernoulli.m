clear
clc

a = 1; % semi-major axis
b = 1; % semi-minor axis
syms t
% Compute x and y using parametric equations
x = a*cos(t);
y = b*sin(t)* cos(t);

% Compute derivatives of x and y with respect to time
dx = diff(x,t);
dy = diff(y,t);
v = sqrt(dx^2 + dy^2);

% Compute angular velocity (rate of change of phi with respect to time)
d2x = diff(dx,t);
d2y = diff(dy,t);
theta = atan2(dy,dx);
w = (dx*d2y-dy*d2x)/(dx^2+dy^2);


%%

num_points = 100;  % number of points
k = 2;
xo = 0;
yo = 0;

% Compute theta values
t = linspace(0, 2*pi*k, num_points);  % linear theta values from 0 to 2*pi;
x = a*cos(t);
y = b*sin(t).*cos(t);

plot(x,y)
title('Lemniscata Trajectory')
figure
subplot(2,1,1)
v = (sin(t).^2 + (cos(t).^2 - sin(t).^2).^2).^(1/2);
plot(t,v)
title('v_{ref}')
subplot(2,1,2)
w = (cos(t).*(cos(t).^2 - sin(t).^2) + 4*cos(t).*sin(t).^2)./(sin(t).^2 + (cos(t).^2 - sin(t).^2).^2);
plot(t,w)
title('w_{ref}')
