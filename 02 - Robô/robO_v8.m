clear;close all;clc

% Inicializações
pvar e1 e2 e4 e5
e = [e1;e2;e4;e5];
ebar = [0;0;0;1];
Ur = [2;2];
ubar = Ur;
% Input-affine Robot tracking error dynamics
f = [Ur(1)*e5+e2*ubar(2); Ur(1)*e4-e1*ubar(2); Ur(2)*e5; 0];
g = [-1 e2; 0 -e1; 0 -e5-1; 0 e4];
% Vector of monomials chosen
Z = e;
% Linear-like Robot tracking error dynamics 
A =  [0 ubar(2) 0 Ur(1);
      -ubar(2) 0 Ur(1) 0;
      0 0 0 Ur(2);
      0 0 0 0];
B = g;
% Dimensions
n = length(e);
m = size(g,2);
nz = length(Z);

% definicao de X = E(Sx)
epsi = 1e-8;
% -- SOS PROGRAM
prog = sosprogram(e);
[prog,Q] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');
[prog,K] = sospolymatrixvar(prog,monomials(e,0),[m,nz]);
[prog,l1] = sossosvar(prog,monomials(e,1:2));
[prog,l2] = sossosvar(prog,monomials(e,1:2));
[prog,s1] = sossosvar(prog,monomials(e,1:2));
[prog,s2] = sossosvar(prog,monomials(e,1:2));
%geq = e4^2+e5^2+2*e5; % restrição de igualdade
geq = e4^2+e5^2-1;
for i=1:nz
    for j=1:n
    M(i,j) = diff(Z(i),e(j)); % matriz M (derivadas parciais de Z(x))
    end
end
%prog = sosineq(prog,l1-1e-8);
%prog = sosineq(prog,l2-1e-8*(e'*e));
S1 = Z'*Q*Z-geq*s1-l1;
%S1 = Q+geq*Seq1*eye(nz)-epsi*eye(nz);
prog = sosineq(prog,S1);
S2 = -Z'*(Q*A'*M'+M*A*Q+K'*B'*M'+M*B*K)*Z-geq*s2-l2;
prog = sosineq(prog,S2);
% Solution
sol = sossolve(prog);
Q = double(sosgetsol(sol,Q));
K = sosgetsol(sol,K);
l1 = sosgetsol(sol,l1);
l2 = sosgetsol(sol,l2);
s1 = sosgetsol(sol,s1);
s2 = sosgetsol(sol,s2);
% Lyapunov and Controllers
V = Z'*inv(Q)*Z;
F = K*inv(Q);
u = F*Z;

% de/dt
f = [Ur(1)*e5; Ur(1)*e4; Ur(2)*e5; -Ur(2)*e4];
g = [-1 e2; 0 -e1; 0 -e5; 0 e4];
global fx gx ux
fx = matlabFunction(p2s(f));
gx = matlabFunction(p2s(g));
ux = matlabFunction(p2s(u));
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
close all
%e0 = [-1.802; 3.133; sin(0.791);cos(0.791)];
xr0 = [1;1;0];
xc0 = [0;0;0];
th = xc0(3);
e0 = xr0-xc0;
R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
e0 = R'*e0;
e0 = [e0(1:2);sin(e0(3));cos(e0(3));xr0;xc0];

opts = odeset('MaxStep',1e-1);
[t, e] = ode45(@closedloop, [0 10], e0, opts);
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

% Plot do sinal de controle
u1 = matlabFunction(p2s(u(1)));
u2 = matlabFunction(p2s(u(2)));
figure
response_U1 = plot(t,u1(e(:,1),e(:,2),e(:,3),e(:,4)-1)+ubar(1));
hold on
response_U2 = plot(t,u2(e(:,1),e(:,2),e(:,3),e(:,4)-1)+ubar(2));
legend('v','w')

function de = closedloop(t,e)
    global fx gx ux gxr gxc
    de = [fx(e(3),e(4))+gx(e(1),e(2),e(3),e(4))*(ux(e(1),e(2),e(3),e(4)-1)+[2;2]); % de/dt
          gxr(e(7))*[2;2]; %dxr/dt
          gxc(e(10))*(ux(e(1),e(2),e(3),e(4)-1)+[2;2])];
end





