clc; clear
close all

T = 100; % samples
ts = 0.5;  % sample time
delta = 0.1;

A = [0 1; 0 0];
B = [0 ; 1];
A_dt=eye(2)+ts*A;   % euler discretization of A
B_dt=ts*B;

n = size(A,2); % colunas de A
m = size(B,2); % colunas de B

DELTA = sqrt(T*delta)*eye(n);

Uo = 2*rand(m,T)-1;
Xo = rand(n,1);
for i=2:T+1
    d = sqrt(delta)*[cos(2*pi*0.4*(ts*(i-1))) ; sin(2*pi*0.4*(ts*(i-1)))];
    Xo(:,i)=A_dt*Xo(:,i-1)+B_dt*Uo(i-1)+d;
end

X1=Xo(:,2:T+1);
Xo=Xo(:,1:T);

subplot(2,1,1)
plot(0:T-1,Uo);
title('Entrada')
subplot(2,1,2)
plot(0:T-1,Xo(1,:),0:T-1,Xo(2,:))
title('Estados')
legend('x1','x2','Location','best')

Ad = [Xo;Uo]*[Xo;Uo]';
Bd = -[Xo;Uo]*X1';
Cd = X1*X1'-DELTA*DELTA';

Y = sdpvar(m,n,'full','real');
P = sdpvar(n,n,'symmetric','real');
CE = [];
CE = CE+(P>=0);
CE = CE+([-P-Cd zeros(n,n) Bd'; zeros(n,n) -P [P;Y]'; Bd [P;Y] -Ad]<=0);

options = sdpsettings('solver','sedumi','verbose',0);
optimize(CE);

Y = value(Y);
P = value(P);

K = Y*inv(P)

for i=2:T+1
    Xo(:,i)=A_dt*Xo(:,i-1)+B_dt*K*Xo(:,i-1);
end

figure
plot(0:T,Xo(1,:),0:T,Xo(2,:))
title('Estados')
legend('x1','x2','Location','best')

figure
a=4;
passo = 0.01;
[x1,x2] = meshgrid(-a:passo:a,-a:passo:a);
%malha aberta
dx1 = A(1,1).*x1+A(1,2).*x2;
dx2 = A(2,1).*x1+A(2,2).*x2;
streamslice(x1,x2,dx1,dx2)
title('Malha aberta'); xlabel('x1'); ylabel('x2');
figure
%malha fechada
dx1 = A(1,1).*x1+A(1,2).*x2;
dx2 = A(2,1).*x1+A(2,2).*x2+K(1).*x1+K(2).*x2;
streamslice(x1,x2,dx1,dx2)
title('Malha fechada'); xlabel('x1'); ylabel('x2');

pvar x1 x2;
x = [x1;x2];
V = x'*P*x

