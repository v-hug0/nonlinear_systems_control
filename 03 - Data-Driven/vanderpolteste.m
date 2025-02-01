pvar x1 x2
vars = [x1;x2];
x = [x1;x2];
e = 1;
% f = [x2;-x1+e*(1-x1^2)*x2]
% g = [1 1];
% h = x1+x2+x1^2+x1*x2+x2^2;
Z = monomials(x,1:3);
A = [0 1 0 0 0 0 0 0 0; -1 e 0 0 0 0 -1 0 0];
W = monomials(x,0);
B = [1;1];

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
d = zeros(n,T);

for i=2:T+1
     d(:,i-1) = sqrt(delta)*[cos(2*pi*0.4*t(i-1));sin(2*pi*0.4*t(i-1))]; 
     Z0(:,i-1) = [X0(1,i-1);X0(2,i-1);X0(1,i-1)^2;X0(1,i-1)*X0(2,i-1);X0(2,i-1)^2;X0(1,i-1)^3;X0(1,i-1)^2*X0(2,i-1);X0(1,i-1)*X0(2,i-1)^2;X0(2,i-1)^3];
     W0(:,i-1) = [1];
     V0(:,i-1) = W0(:,i-1)*u(i-1);
     X1(:,i-1) = A*Z0(:,i-1)+B*V0(:,i-1)+d(:,i-1);
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









