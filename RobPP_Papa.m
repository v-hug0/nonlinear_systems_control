clear;close all;clc

%%% Initializations
pvar e1 e2
e=[e1;e2];

% Polynomial recasted robot tracking model from eq.(16).
g = [-1 e2; 0 -e1];

% Input-affine Robot tracking error dynamics
f = [0;0];
g = [-1 e2; 0 -e1];
% Vector of monomials chosen
Z = e;
nz = length(Z);
% Linear-like Robot tracking error dynamics 
% Dimensions
n = length(e);
m = size(g,2);

A = zeros(n,nz);
B = g;

% definicao de X = E(Sx)
epsi = 1e-8;
% -- SOS PROGRAM
prog = sosprogram(e);
[prog,Q] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');
[prog,K] = sospolymatrixvar(prog,monomials(e,0:3),[m,nz]);
[prog,l1] = sospolyvar(prog,monomials(e,2:4));
[prog,l2] = sospolyvar(prog,monomials(e,2:6));
for i=1:nz
    for j=1:n
    M(i,j) = diff(Z(i),e(j)); % matriz M (derivadas parciais de Z(x))
    end
end
S1 = Q-l1;
prog = sosineq(prog,S1);
S2 = -(Q*A'*M'+M*A*Q+K'*B'*M'+M*B*K)-l2;
prog = sosineq(prog,S2);
% Solution
sol = sossolve(prog);
Q = double(sosgetsol(sol,Q));
K = sosgetsol(sol,K);
% Lyapunov and Controllers
V = Z'*inv(Q)*Z;
F = K*inv(Q);
u = F*Z;

% de/dt
g = [-1 e2; 0 -e1];
global gx ux
gx = matlabFunction(p2s(g));
ux = matlabFunction(p2s(u));
syms xc3
% dxr/dt
gc = [cos(xc3) 0; sin(xc3) 0; 0 1];
%gc = [e5 0; e4 0; 0 e5; 0 -e4];
global gxc
gxc = matlabFunction(gc);

close all
xr0 = [-1;1;0];
xc0 = [0;0;0];
th = xc0(3);
e0 = xr0-xc0;
R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
e0 = R'*e0;
e0 = [e0(1:2);xc0];

opts = odeset('MaxStep',1e-1);
[t, e] = ode45(@closedloop, [0 10], e0, opts);
response_CL = plot(t,e(:,1),t,e(:,2));
legend('e1','e2','e4','e5')

% Plot das trajetórias
figure
real_trajectory = plot(e(:,3),e(:,4));
hold on
real0_trajectory = plot(e(end,3),e(end,4),'o');

% Plot do sinal de controle
u1 = matlabFunction(p2s(u(1)));
u2 = matlabFunction(p2s(u(2)));
figure
response_U1 = plot(t,u1(e(:,1),e(:,2)));
hold on
response_U2 = plot(t,u2(e(:,1),e(:,2)));
legend('v','w')

function de = closedloop(t,e)
    global gx ux gxr gxc
    de = [gx(e(1),e(2))*ux(e(1),e(2)); % de/dt
          gxc(e(5))*ux(e(1),e(2))];
end





