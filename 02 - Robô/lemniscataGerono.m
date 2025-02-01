clear
clc

a = 0.4;
syms t
% Compute x and y using parametric equations
x = a*cos(t)/(1+sin(t)^2);
y = a*sin(t)*cos(t)/(1+sin(t)^2);

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
x = a*cos(t)./(1+sin(t).^2) + xo;
y = a*sin(t).*cos(t)./(1+sin(t).^2) + yo;

plot(x,y)
title('Lemniscata Trajectory')
figure
subplot(2,1,1)
v = (((2.*sin(t))./(5.*(sin(t).^2 + 1)) + (4.*cos(t).^2.*sin(t))./(5.*(sin(t).^2 + 1).^2)).^2 + ((2.*sin(t).^2)./(5.*(sin(t).^2 + 1)) - (2.*cos(t).^2)./(5.*(sin(t).^2 + 1)) + (4.*cos(t).^2.*sin(t).^2)./(5.*(sin(t).^2 + 1).^2)).^2).^(1./2);
plot(t,v)
title('v_{ref}')
subplot(2,1,2)
w = -(((2.*sin(t).^2)./(5.*(sin(t).^2 + 1)) - (2.*cos(t).^2)./(5.*(sin(t).^2 + 1)) + (4.*cos(t).^2.*sin(t).^2)./(5.*(sin(t).^2 + 1).^2)).*((2.*cos(t))./(5.*(sin(t).^2 + 1)) + (4.*cos(t).^3)./(5.*(sin(t).^2 + 1).^2) - (12.*cos(t).*sin(t).^2)./(5.*(sin(t).^2 + 1).^2) - (16.*cos(t).^3.*sin(t).^2)./(5.*(sin(t).^2 + 1).^3)) - ((2.*sin(t))./(5.*(sin(t).^2 + 1)) + (4.*cos(t).^2.*sin(t))./(5.*(sin(t).^2 + 1).^2)).*((8.*cos(t).*sin(t))./(5.*(sin(t).^2 + 1)) - (12.*cos(t).*sin(t).^3)./(5.*(sin(t).^2 + 1).^2) + (12.*cos(t).^3.*sin(t))./(5.*(sin(t).^2 + 1).^2) - (16.*cos(t).^3.*sin(t).^3)./(5.*(sin(t).^2 + 1).^3)))./(((2.*sin(t))./(5.*(sin(t).^2 + 1)) + (4.*cos(t).^2.*sin(t))./(5.*(sin(t).^2 + 1).^2)).^2 + ((2.*sin(t).^2)./(5.*(sin(t).^2 + 1)) - (2.*cos(t).^2)./(5.*(sin(t).^2 + 1)) + (4.*cos(t).^2.*sin(t).^2)./(5.*(sin(t).^2 + 1).^2)).^2);
plot(t,w)
title('w_{ref}')
