clear;close all;clc

pvar x1 x2 x3 x4
x = [x1;x2;x3;x4];

m = 0.231;
L = 0.32;
Mt = 1.142;
Jm = 0.03153;
g = 9.81;
Mbar=[Mt  m*L; m*L Jm];
c1 = 1.329;
c2 = 5.561;

f = [x3;x4;inv(Mbar)*[m*L*x2*x4^2-c2*x3;m*g*L*x2]];
g = [0;0; inv(Mbar)*[c1;0]];
gw = [0;0;inv(Mbar)*[0;1]];

%%

Z = e;

n = length(e);
m = size(g,2);
nz = length(Z);

A =  [0 0 1 0 0 0 0 0;
      0 0 0 1 0 0 0 0;
      0 -1.755 -5.7407

B = g;

% definition of X = E(Sx)
Sx = eye(nz)*0.1;
epsi = 1e-2;

% -- SOS PROGRAM
prog = sosprogram(e);
[prog,Q] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');
[prog,K] = sospolymatrixvar(prog,monomials(e,0:3),[m,nz]);
[prog,T] = sospolymatrixvar(prog,monomials(e,0:3),[m,nz]);
gx = 1-Z'*Sx*Z; %elipsoid of X
geq = e4^2+e5^2-1; %equality constraint
for i=1:nz
    for j=1:n
    M(i,j) = diff(Z(i),e(j));
    end
end
% -- Theorem 1 - Eq.(11)
S3 = Z'*(Sx-Q)*Z;
prog = sosineq(prog,S3);
% -- Theorem 1 - Eq.(12)
% vertΩ set of θ
th1 = [0;0];
th2 = [0;1];
th3 = [1;0];
th4 = [1;1];
vert_omega = [th1 th2 th3 th4];
n_vert = 2^m; % n° vertices of omega
for i=1:n_vert
    [prog,Sr] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');
    [prog,Seq] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');
    TH = eye(m)*vert_omega(i);
    F1 = -(M*A*Q+M*B*(TH*K+(eye(m)-TH)*T) + Q*A'*M' + (K'*TH'+T'*(eye(m)-TH'))*B'*M') - Sr*gx - Seq*geq;
    S1 = F1 - epsi*eye(nz);
    prog = sosineq(prog,S1);
end
% Theorem 1 - Eq.(13)
for j=1:m
    [prog,Sh] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');
    [prog,Seq] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');
    F2 = Q-Sh*gx-Seq*geq;
    S2 = [1 T(j); T(j)' F2(j)];
    prog = sosineq(prog,S2);
end
obj = trace(Q);
prog = sossetobj(prog,-obj);

% Solution
sol = sossolve(prog);
Q = double(sosgetsol(sol,Q));
K = sosgetsol(sol,K);
T = sosgetsol(sol,T);

% Lyapunov and Controllers
V = Z'*inv(Q)*Z;
F = K*inv(Q);
H = T*inv(Q);
u = F*Z;
v = H*Z;

%% PLOT

a=1;
passo = 0.01;
[e1,e2] = meshgrid(-a:passo:a,-a:passo:a);

dx1 = -x1+x1.*x2+x2.*K.*(x2+x2.^2+x2.^3);
dx2 = x1+2.*x2+x1.^2+x1.^2.*x2+K.*(x2+x2.^2+x2.^3);
streamslice(x1,x2,dx1,dx2,2)
title('Malha fechada'); xlabel('x1'); ylabel('x2');
hold on;

elipV = matlabFunction(p2s(1-Z'*Sx*Z));
zhandle = fimplicit(elipV)
zhandle.LineWidth = 1.5;

elipL = matlabFunction(p2s(1-V));
zhandle = fimplicit(elip)
zhandle.LineWidth = 1.5;
zhandle.Color = "b";

%%
