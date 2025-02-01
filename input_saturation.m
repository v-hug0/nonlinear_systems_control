

pvar x1 x2 v1 v2
var = [x1;x2;v1;v2];
x = [x1;x2];
v = [v1;x2];

g = 17.76 - 103.79*x1 + 229.62*x1^2 - 226.31*x1^3 + 83.72*x1^4;
A = [-0.5*g 0.5; -0.2 -0.3];
B = [0;0.2];

u_nom = 1.2;
x1bar = 1;
x2bar=(17.76 - 103.79*x1bar + 229.62*x1bar^2 - 226.31*x1bar^3 + 83.72*x1bar^4);
xbar=[x1bar;x2bar];
ubar=(-0.2*x1bar-0.3*x2bar)/0.2;


%% malha aberta

syms x1 x2;
a = 2;
[x1,x2] = meshgrid(-a:0.01:a,-a:0.01:a);
dx1 = -0.5*(17.76 - 103.79.*x1 + 229.62.*x1.^2 - 226.31.*x1.^3 + 83.72.*x1.^4).*x1 + A(1,2).*x2;
dx2 = A(2,1).*x1 + A(2,2).*x2 + 0.2*u_nom;
figure
streamslice(x1,x2,dx1,dx2,1)
title('Diagrama de fases - malha aberta')
xlabel('x1'),ylabel('x2')


%% malha fechada

syms x1 x2;
a = 4;
[x1,x2] = meshgrid(-a:0.01:a,-a:0.01:a);
Kx = K(1).*(x1-x1bar) + K(2).*(x2-x2bar);
dx1 = -0.5*(17.76 - 103.79.*(x1-x1bar) + 229.62.*(x1-x1bar).^2 - 226.31.*(x1-x1bar).^3 + 83.72.*(x1-x1bar).^4).*(x1-x1bar) + A(1,2).*(x2-x2bar) + B(1)*Kx;
dx2 = A(2,1).*(x1-x1bar) + A(2,2).*(x2-x2bar) + B(2)*Kx;
figure
streamslice(x1,x2,dx1,dx2,1)
title('Diagrama de fases - malha fechada')
xlabel('x1'),ylabel('x2')
hold on
plot(x1bar,x2bar,'o')

%%


n = size(A,1);
nz = length(Z);

beta1 = 1e-4;

prog = sosprogram(var);
[prog,P] = sospolymatrixvar(x,monomials(x,0),[nz,nz]);
[prog,K] = sospolymatrixvar(x,monomials(x,4),[n,nz]);
[prog,beta1] = sospolyvar(x,monomials(x,2:4));
s1 = v'*(P-beta1*eye(nz))*v;
s2 = -v'*(P*A'*M'+M*A*P+K'*B'*M'+M*B*K+beta2*eye(nz))*v;
prog = sosineq(prog,s1);
prog = sosineq(prog,s2);
sol = sossolve(prog);
P = sosgetsol(P);
K = sosgetsol(K);

V = Z'*inv(P)*Z;
F = K*inv(P);
u = F*Z;

