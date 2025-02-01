  clc; clear; close all

% parametros do circuito

L = 2e-3;
Cap = 220e-6;
R = 3;
Vin = 12;

x1bar = 0.8333;
x2bar = 2.5;
ubar = 0.2083;

%

pvar x1 x2
vars = [x1;x2];
A = [0 -1/L; 1/Cap -1/(R*Cap)];
B = [Vin/L; 0];
C = [1 0; 0 1]; %%% como definir?
D = [0 0; 0 0]; %%% como definir?

n = size(A,2); % tamanho de x
m = size(B,2); % tamanho de u

T = 50; % samples
delta = 1e-4;
DELTA = T*delta*eye(n);

tsim = 0.01;
tsample = tsim/T;
t = 0:tsample:tsim;
w = 2*pi*(3/tsim);
U0 = 0.25*sin(w*t+pi/4)+0.5;

% gerando os dados a partir do modelo

X0 = zeros(n,T+1);
X0(:,1) = 1;
d = zeros(n,T);

for i=2:T+1
    d(:,i-1) = sqrt(delta)*[cos(2*pi*1e-3*t(i-1)); sin(2*pi*1e-3*t(i-1))];
    X0(:,i) = X0(:,i-1) + tsample*(A*X0(:,i-1)+B*U0(i-1)+d(:,i-1));
end

X1 = (X0(:,2:T+1) - X0(:,1:T))/tsample;

% translação

X0 = X0(:,1:T)+[x1bar;x2bar];
U0 = U0(1:T)+ubar;

% grafico dos dados gerados

t = t(1:T);
figure
subplot(2,1,1)
plot(t,U0)
title('Entrada')

subplot(2,1,2)
plot(t,X0(1,:),t,X0(2,:))
legend('x1','x2','Location','best')
title('Malha aberta')


% matrizes de dados

Ad = [X0;U0]*[X0;U0]';
Bd = -[X0;U0]*X1';
Cd = X1*X1'-DELTA*DELTA';
zeta = -inv(Ad)*Bd;
Q = Bd'*inv(Ad)*Bd-Cd;

% declaração das variaveis LMI

prog = sosprogram(vars);
[prog,Y] = sospolymatrixvar(prog,monomials(vars,0),[m,n]);
[prog,P] = sospolymatrixvar(prog,monomials(vars,0),[n,n],'symmetric');

% condicoes de estabilidade

prog = sosineq(prog,P);
cond = [-Cd Bd'-[P;Y]'; Bd-[P;Y] -Ad];
prog = sosineq(prog,-cond);

% resolvendo

sol = sossolve(prog);
Y = double(sosgetsol(sol,Y));
P = double(sosgetsol(sol,P));

% controlador

K = Y*inv(P);

% funcao de Lyapunov

pvar x1 x2
x = [x1;x2];
V = x'*P*x

% resposta do sistema em MF

buck_mf = ss(A,B*K,C,D);
[X0_mf,tOut] = initial(buck_mf,X0(:,1),t);
plot(t,X0_mf'+[x1bar;x2bar])
legend('x1','x2','Location','best')
title('K (projetado)')

% diagrama de fases - malha aberta

syms x1 x2;
a = 4;
[x1,x2] = meshgrid(-a:0.01:a,-a:0.01:a);
dx1 = A(1,1).*x1 + A(1,2).*x2;
dx2 = A(2,1).*x1 + A(2,2).*x2;

figure
streamslice(x1,x2,dx1,dx2,1)
title('Diagrama de fases - malha aberta')
xlabel('x1'),ylabel('x2')

% diagrama de fases - malha fechada

syms x1 x2;
a = 4;
[x1,x2] = meshgrid(-a:0.01:a,-a:0.01:a);
Kx = K(1).*(x1-x1bar) + K(2).*(x2-x2bar) + ubar;
dx1 = A(1,1).*(x1-x1bar) + A(1,2).*(x2-x2bar) + B(1)*Kx;
dx2 = A(2,1).*(x1-x1bar) + A(2,2).*(x2-x2bar) + B(2)*Kx;

figure
streamslice(x1,x2,dx1,dx2,1)
title('Diagrama de fases - malha fechada')
xlabel('x1'),ylabel('x2')
hold on
plot(x1bar,x2bar,'o')


