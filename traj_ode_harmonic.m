clear
clc

a = 0.4; % 400 mm/s;
b = 0.6; % 600 mm
wi = 2*pi/20;
%w = 2*pi/5;
syms t
x = a*t;

y = b*sin(wi*t);

dx = diff(x,t);
dy = diff(y,t);

ddx = diff(dx,t);
ddy = diff(dy,t);
global v
v = matlabFunction(sqrt(dx^2+dy^2));
global w
w = matlabFunction((dx*ddy-dy*ddx)/(dx^2+dy^2));

%%

pvar xr4 xr5
% dxr/dt
gr = [xr5 0; xr4 0; 0 xr5; 0 -xr4];
global gxr
gxr = matlabFunction(p2s(gr));

xr0 = [0;0;0];
xr0 = [xr0(1:2);sin(xr0(3));cos(xr0(3))];
opts = odeset('MaxStep',1e-1);
tsim = [0 30];
% Plot da trajetória de referência do robô
[t, xr] = ode45(@ref_rob, tsim, xr0, opts);
figure
ref_trajectory = plot(xr(:,1),xr(:,2));
hold on
refF_trajectory = plot(xr(1,1),xr(1,2),'o');
figure
plot(t,subs(v),t,subs(w))
legend('v','w')

function dxr = ref_rob(t,xr)
    global gxr v w
    dxr = gxr(xr(3),xr(4))*[v(t);w(t)];
end







