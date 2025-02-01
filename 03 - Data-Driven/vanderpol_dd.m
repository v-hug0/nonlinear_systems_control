
clc; clear; close all

pvar x1 x2
vars = [x1;x2];
x = [x1;x2];
e = 1;

A = [0 1 0 0 0 0 0 0 0; -1 e 0 0 0 0 -1 0 0];
Z = monomials(x,1:3);
B = [1;1];
W = monomials(x,0);

n = 2; % tamanho de x
m = 1; % tamanho de u

N = length(Z);
M = length(W);

T = 100;
delta = 1e-5;
DELTA = T*delta*eye(n);
tsim = 1;
tsample = tsim/T;
t = 0:tsample:tsim;
w = 2*pi*(3/tsim);
u = 0.5*(sin(w*t)+1);

% gerando dados

X0 = zeros(n,T+1);
X0(:,1) = 0.1;
Z0 = zeros(N,T);
W0 = zeros(M,T);
V0 = zeros(M,T);
X1 = zeros(n,T);
d = zeros(n,T+1);

for i=2:T+1
    d(:,i-1) = sqrt(delta)*[cos(2*pi*0.4*t(i-1));sin(2*pi*0.4*t(i-1))];
    Z0(:,i-1) = [X0(1,i-1); X0(2,i-1); X0(1,i-1)^2; X0(1,i-1)*X0(2,i-1);X0(2,i-1)^2; X0(1,i-1)^3; X0(1,i-1)^2*X0(2,i-1);X0(1,i-1)*X0(2,i-1)^2; X0(2,i-1)^3];
    W0(:,i-1) = [1];
    V0(:,i-1) = W0(:,i-1)*u(i-1);
    X1(:,i-1) = A*Z0(:,i-1)+B*V0(:,i-1) +d(:,i-1);
    X0(:,i) = X0(:,i-1)+tsample*X1(:,i-1);
end

X0 = X0(:,1:T);
V0 = V0(:,1:T);
t = t(1:T);
u = u(1:T);
 
% grafico dos dados gerados
 
figure
subplot(2,1,1)
plot(t,u)
title('Entrada')

subplot(2,1,2)
plot(t,X0(1,:),t,X0(2,:))
legend('x1','x2','Location','best')
title('Malha aberta')

% matrizes de dados
 
Ad = [Z0;V0]*[Z0;V0]';
Bd = -[Z0;V0]*X1';
Cd = X1*X1'-DELTA*DELTA';
 
zeta = -inv(Ad)*Bd;
Q = Bd'*inv(Ad)*Bd-Cd;

% parametros de projeto

c = 1;
epsi = 1e-5;
l0 = 1e0*(x1^2+x2^2);

% V inicial 

V = x'*x;

i = 0;
i_max = 15;
%% 

while i<i_max
    % solving with SOSTOOLS (V fixed)
    prog = sosprogram(vars);
    
    [prog,k] = sospolyvar(prog,[monomials(x,1:2)],'wscoeff');
    [prog,lambda] = sospolyvar(prog,[monomials(x,0:2)],'wscoeff');
    [prog,Tx] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    [prog,l1] = sospolyvar(prog,[monomials(x,2)],'wscoeff');
    [prog,s1] = sospolyvar(prog,[monomials(x,2)],'wscoeff');
    [prog,s2] = sospolymatrixvar(prog,[monomials(x,0)],[m+N+M+n,m+N+M+n]);
    
    % SOS constraints
    
    
    prog = sosineq(prog,V-l1);
    
    gradV = [diff(V,x1) diff(V,x2)];
    estabilidade = [Tx+gradV*zeta'*[Z;W*k] [Z;W*k]'*Ad^(-1/2) gradV*lambda*Q^(1/2);
                    Ad^(-1/2)*[Z;W*k] -lambda*eye(N+M) zeros(N+M,n);
                    lambda*Q^(1/2)*gradV' zeros(n,N+M) -4*lambda*eye(n)];
    prog = sosineq(prog,-estabilidade);
    
    prog = sosineq(prog,lambda-epsi);
    
    % first solution (k and lambda)
    sol = sossolve(prog);
    k = sosgetsol(sol,k);
    lambda = sosgetsol(sol,lambda);
    
    % solving with SOSTOOLS (k and lamda fixed)
    prog = sosprogram(vars);
    
    [prog,V] = sospolyvar(prog,[monomials(x,2)],'wscoeff');
    [prog,Tx] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    [prog,l1] = sospolyvar(prog,[monomials(x,2)],'wscoeff');
    [prog,s1] = sospolyvar(prog,[monomials(x,2)],'wscoeff');
    [prog,s2] = sospolymatrixvar(prog,[monomials(x,0)],[m+N+M+n,m+N+M+n]);
    
    % SOS constraints
    prog = sosineq(prog,V-l1);
    
    gradV = [diff(V,x1) diff(V,x2)];
    estabilidade = [Tx+gradV*zeta'*[Z;W*k] [Z;W*k]'*Ad^(-1/2) gradV*lambda*Q^(1/2);
                    Ad^(-1/2)*[Z;W*k] -lambda*eye(N+M) zeros(N+M,n);
                    lambda*Q^(1/2)*gradV' zeros(n,N+M) -4*lambda*eye(n)];
    prog = sosineq(prog,-estabilidade);
    prog = sosineq(prog,lambda-epsi);
    
    % first solution (V)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    i = i+1
end

%% 
Kc = full(k.coefficient);

for i=2:T+1
    Z0(:,i-1) = [X0(1,i-1); X0(2,i-1); X0(1,i-1)^2; X0(1,i-1)*X0(2,i-1);
                 X0(2,i-1)^2; X0(1,i-1)^3; X0(1,i-1)^2*X0(2,i-1);
                 X0(1,i-1)*X0(2,i-1)^2; X0(2,i-1)^3];
    W0(:,i-1) = [1];

    K = Kc(1)*X0(1,i-1)^2 + Kc(2)*X0(1,i-1)*X0(2,i-1) + Kc(3)*X0(2,i-1)^2 +...
    Kc(4)*X0(1,i-1) + Kc(5)*X0(2,i-1);
       
    X0(:,i) = X0(:,i-1)+tsample*(A*Z0(:,i-1)+B*W0(:,i-1)*K+d(:,i-1));
end

figure
plot(t,X0(1,1:T),t,X0(2,1:T))
legend('x1','x2','Location','best')
title('Malha fechada')

 % diagrama de fases - malha aberta

syms x1 x2
a = 4;
[x1,x2] = meshgrid(-a:0.01:a,-a:0.01:a);
dx1 = x2;
dx2 = -x1 + e.*(1-x1.^2).*x2;
figure
streamslice(x1,x2,dx1,dx2,2)
title('Diagrama de fases - malha aberta')
xlabel('x1'),ylabel('x2')


% diagrama de fases - malha fechada

syms x1 x2
a=0.5;
[x1,x2] = meshgrid(-a:0.01:a,-a:0.01:a);
K = Kc(1).*x1.^2 + Kc(2).*x1.*x2 + Kc(3).*x2.^2 + Kc(4).*x1 + Kc(5).*x2;
dx1 = x2 + K;
dx2 = -x1 + e.*(1-x1.^2).*x2 + K;

figure
streamslice(x1,x2,dx1,dx2,4)
title('Diagrama de fases - malha fechada')
xlabel('x1'),ylabel('x2')
hold on
plot(0,0,'o')






