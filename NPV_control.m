clear;close all;clc
%%
% Inicializações - CIRCULO
r = 1;
w = 2;
syms t
x = r*cos(w*t);
y = r*sin(w*t);
dx = diff(x,t);
dy = diff(y,t);
ddx = diff(dx,t);
ddy = diff(dy,t);

global v w dv dw

v = sqrt(dx^2+dy^2);
w = (dx*ddy-dy*ddx)/(dx^2+dy^2);
dv = diff(v,t);
dw = diff(w,t);
v = matlabFunction(v);
w = matlabFunction(w);
dv = matlabFunction(dv);
dw = matlabFunction(dw);

%%

% Inicializações - LEMNISCATA
l = 2.5; % 2500 mm
w = pi/10;
syms t
x = l*sin(2*w*t);
y = l*cos(w*t);

dx = diff(x,t);
dy = diff(y,t);

ddx = diff(dx,t);
ddy = diff(dy,t);

global v w dv dw

v = sqrt(dx^2+dy^2);
w = (dx*ddy-dy*ddx)/(dx^2+dy^2);
dv = diff(v,t);
dw = diff(w,t);
v = matlabFunction(v);
w = matlabFunction(w);
dv = matlabFunction(dv);
dw = matlabFunction(dw);


%%

% Inicializações - HARMONICA
a = 0.4; % 400 mm/s;
b = 0.6; % 600 mm
wi = pi/5;
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

pvar e1 e2 e4 e5
e = [e1;e2;e4;e5];
ebar = [0;0;0;1];
global ubar
t = 0;
pvar vr wr dvr dwr
sig = [vr;wr];
dsig = [dvr;dwr];
ubar = sig;
% Input-affine Robot tracking error dynamics
f = [sig(1)*e5+e2*ubar(2); sig(1)*e4-e1*ubar(2); sig(2)*e5; 0];
g = [-1 e2; 0 -e1; 0 -e5-1; 0 e4];
% Vector of monomials chosen
Z = e;
% Linear-like Robot tracking error dynamics 
A =  [0 ubar(2) 0 sig(1);
      -ubar(2) 0 sig(1) 0;
      0 0 0 sig(2);
      0 0 0 0];
B = g;
% Dimensions
n = length(e);
m = size(g,2);
nz = length(Z);

% definicao de X = E(Sx)
epsi = 1e-8;
% -- SOS PROGRAM
prog = sosprogram([e;sig;dsig]);
[prog,Q] = sospolymatrixvar(prog,monomials(sig,0:1),[nz,nz],'symmetric');
[prog,K] = sospolymatrixvar(prog,monomials([e;sig;dsig],0:1),[m,nz]);
[prog,l1] = sospolyvar(prog,monomials(e,2:4));
[prog,l2] = sospolyvar(prog,monomials(e,2:4));
[prog,s1] = sospolyvar(prog,monomials(e,1:2));
[prog,s2] = sospolyvar(prog,monomials(e,1:2));
%geq = e4^2+e5^2+2*e5; % restrição de igualdade
geq = e4^2+e5^2-1;
for i=1:nz
    for j=1:n
    M(i,j) = diff(Z(i),e(j)); % matriz M (derivadas parciais de Z(x))
    end
end
S1 = Q-geq*s1-l1*eye(nz);
%S1 = Q+geq*Seq1*eye(nz)-epsi*eye(nz);
prog = sosineq(prog,S1);
S2 = -(Q*A'*M'+M*A*Q+K'*B'*M'+M*B*K-(diff(Q,sig(1))*dsig(1)+diff(Q,sig(2))*dsig(2)))-geq*s1-l2;
prog = sosineq(prog,S2);
% Solution
sol = sossolve(prog);
Q = sosgetsol(sol,Q);
K = sosgetsol(sol,K);

% Lyapunov and Controller
global Q1 K1
Q1 = matlabFunction(p2s(Q));
K1 = matlabFunction(p2s(K));
%%
% de/dt
f = [sig(1)*e5; sig(1)*e4; sig(2)*e5; -sig(2)*e4];
g = [-1 e2; 0 -e1; 0 -e5; 0 e4];
global fx gx ux
fx = matlabFunction(p2s(f));
gx = matlabFunction(p2s(g));

syms xr3
% dxr/dt
gr = [cos(xr3) 0; sin(xr3) 0; 0 1];
%gr = [e5 0; e4 0; 0 e5; 0 -e4];
global gxr
gxr = matlabFunction(gr);
syms xc3
% dxr/dt
gc = [cos(xc3) 0; sin(xc3) 0; 0 1];
%gc = [e5 0; e4 0; 0 e5; 0 -e4];
global gxc
gxc = matlabFunction(gc);

%%
%e0 = [-1.802; 3.133; sin(0.791);cos(0.791)];
close all
xr0 = [0;0;0];
xc0 = [-1.6;0.3;0];
%xc0 = [1;2;0];
deg0 = 180*xc0(3)/pi
e0 = xr0-xc0;
th=xc0(3);
R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1]; % rotation matrix
e0 = R'*e0;
e0 = [e0(1:2);sin(e0(3));cos(e0(3));xr0;xc0];
%xc0 = [xc0(1:2);sin(xc0(3));cos(xc0(3))];
%xr0 = [xr0(1:2);sin(xr0(3));cos(xr0(3))];
opts = odeset('MaxStep',1e-1);
[t, e] = ode45(@closedloop, [0 20], e0, opts);
response_CL = plot(t,e(:,1),t,e(:,2),t,e(:,3),t,e(:,4));
legend('e1','e2','e4','e5')

% Plot das trajetórias
figure
ref_trajectory = plot(e(:,5),e(:,6));
hold on
ref0_trajectory = plot(e(end,5),e(end,6),'o');
hold on
real_trajectory = plot(e(:,8),e(:,9));
hold on
real0_trajectory = plot(e(end,8),e(end,9),'o');

%%
figure
response_U1 = plot(t,u1(e(:,1),e(:,2),e(:,3),e(:,4)-1)+ubar(1));
hold on
response_U2 = plot(t,u2(e(:,1),e(:,2),e(:,3),e(:,4)-1)+ubar(2));
legend('v','w')
%%
function de = closedloop(t,e)
    global fx gx gxr gxc v w dv dw Q1 K1
    u_ref = [v(t);w(t)];

    K2 = K1(dv(t),dw(t),e(1),e(2),e(3),e(4)-1,v(t),w(t));
    Q1_ = Q1(v(t),w(t));
    Q1_ = inv(Q1_);
    u = K2*inv(Q1_)*[e(1);e(2);e(3);e(4)-1];

    de = [fx(e(3),e(4),v(t),w(t))+gx(e(1),e(2),e(3),e(4))*(u+u_ref); % de/dt
          gxr(e(7))*u_ref; %dxr/dt
          gxc(e(10))*(u+u_ref)];
end

%%

