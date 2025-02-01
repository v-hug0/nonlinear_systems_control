clear;close all;clc

% Inicializações
pvar e1 e2 e4 e5 v1 v2 v3 v4
e = [e1;e2;e4;e5];
v = [v1;v2;v3;v4];
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
var = [e;v];
prog = sosprogram(e);
[prog,Q] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');
[prog,K] = sospolymatrixvar(prog,monomials(e,0:1),[m,nz]);
[prog,l1] = sospolyvar(prog,monomials(e,2:4));
[prog,l2] = sospolyvar(prog,monomials(e,2:4));
[prog,Seq1] = sospolyvar(prog,monomials(e,1:2));
[prog,Seq2] = sospolyvar(prog,monomials(e,1:2));
geq = e4^2+e5^2-1; % restrição de igualdade
for i=1:nz
    for j=1:n
    M(i,j) = diff(Z(i),e(j)); % matriz M (derivadas parciais de Z(x))
    end
end
S1 = Q+geq*Seq1*eye(nz)-l1*eye(nz);
%S1 = Q+geq*Seq1*eye(nz)-epsi*eye(nz);
prog = sosineq(prog,S1);
S2 = -(Q*A'*M'+M*A*Q+K'*B'*M'+M*B*K)+geq*Seq2+l2;
prog = sosineq(prog,S2);
% Solution
sol = sossolve(prog);
Q = double(sosgetsol(sol,Q));
K = sosgetsol(sol,K);
% Lyapunov and Controllers
V = Z'*inv(Q)*Z;
F = K*inv(Q);
u = F*Z;


%% Plot no tempo dos erros

f = [Ur(1)*e5; Ur(1)*e4; Ur(2)*e5; -Ur(2)*e4];
g = [-1 e2; 0 -e1; 0 -e5; 0 e4];
global fx gx ux
fx = matlabFunction(p2s(f));
gx = matlabFunction(p2s(g));
ux = matlabFunction(p2s(u));

e0 = [1; 0; 0; 0];
[t, e] = ode45(@closedloop, [0 15], e0);
response_CL = plot(t,e(:,1),t,e(:,2),t,e(:,3),t,e(:,4));
legend('e1','e2','e4','e5')


function de = closedloop(t,e)
    global fx gx ux
    de = fx(e(3),e(4))+gx(e(1),e(2),e(3),e(4))*(ux(e(1),e(2),e(3),e(4)-1)+[2;2]);
end






%%
% plot
f = [Ur(1)*e5; Ur(1)*e4; Ur(2)*e5; -Ur(2)*e4];
g = [-1 e2; 0 -e1; 0 -e5; 0 e4];
f1 = matlabFunction(p2s(f(1)));
f2 = matlabFunction(p2s(f(2)));
g11 = g(1,1);
g12 = matlabFunction(p2s(g(1,2)));
g21 = g(2,1);
g22 = matlabFunction(p2s(g(2,2)));
u1 = matlabFunction(p2s(u(1)));
u2 = matlabFunction(p2s(u(2)));


%% Elipsoide (não deu certo, vamos tentar plotar no tempo)
a=4;
passo = 0.01;
[e1,e2] = meshgrid(-a:passo:a,-a:passo:a);
ubar=[2;2];
de1 = f1(1)+g11*(u1(e1,e2,0,0)+ubar(1))+g12(e2)*(u2(e1,e2,0,0)+ubar(2));
de2 = f2(0)+g21*(u1(e1,e2,0,0)+ubar(1))+g22(e1)*(u2(e1,e2,0,0)+ubar(2));
streamslice(e1,e2,de1,de2,2)
title('Malha fechada'); xlabel('e1 (x)'); ylabel('e2 (y)');
hold on;
syms e1 e2
elipV = matlabFunction(p2s(1-V));
elipVfix=elipV(e1,e2,0,0);
zhandle = fimplicit(elipVfix)




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           REFERENCE TRAJECTORY GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.01;          % sampling time

Tsim = 20;   %90      % tempo de simulacao
h = dt;           % integrating time


xr = [1;1;0];
Xr = xr;         %Euler


for k = 1:round(Tsim/h)  % 1 a tempo de simulacao/ tempo de integracao   
    vr(k) = Ur(1);
    wr(k) = Ur(2);
    Xr = [Xr xr+h*([vr(k)*cos(xr(3));vr(k)*sin(xr(3));wr(k)])];  %Euler
    xr = Xr(:,k+1);
end
Ur = [vr;wr];
hold on
plot(Xr(1,:),Xr(2,:),'')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          ROBOT REFERENCE TRACKING SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [1;0;0];
x = xo;
X = x;
v = [];
w = [];


u1 = matlabFunction(p2s(u(1)));
u2 = matlabFunction(p2s(u(2)));

for k = 1:(round(Tsim/h))    % 1 a tempo de simulacao/tempo de integracao  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(1);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);
    E5n = E5-1;

    v(k) = u1(E1,E2,E4,E5n) + ubar(1);
    w(k) = u2(E1,E2,E4,E5n) + ubar(2);

    X = [X x+dt*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
end

Error = Xr-X;
Error4 = sin(Error(3,:));
Error5 = cos(Error(3,:));

plot(X(1,:),X(2,:),'r')
legend('Reference Robot','Real Robot')
xlabel('x (m)')
ylabel('y (m)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem de v e w 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
time = dt*(0:Tsim/h-1);
subplot(2,1,1)
plot(time,v,'b')
xlabel('Time (s)')
ylabel('U1 (m/s)')
title('v')
a
hold on
subplot(2,1,2)
plot (time,w,'r')
xlabel('Time (s)')
ylabel('U2 (rad/s)')
title('ω')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem do erro x y e theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(time,Error(1,1:end-1),time,Error(2,1:end-1),time,Error(3,1:end-1),time,Error4(1:end-1),time,Error5(1:end-1));
legend('e1','e2','e3','e4','e5')



