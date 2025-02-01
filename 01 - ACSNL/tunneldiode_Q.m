%% Example 4.1 - Tunnel Diode
clear; close all; clc;

pvar x1 x2
vars = [x1;x2];
x = [x1;x2];
g = 17.76 - 103.79*x1 + 229.62*x1^2 - 226.31*x1^3 + 83.72*x1^4;
A = [-0.5*g 0.5; -0.2 -0.3];
B = [0;0.2];

x1bar = 1;
x2bar=(17.76 - 103.79*x1bar + 229.62*x1bar^2 - 226.31*x1bar^3 + 83.72*x1bar^4);
xbar=[x1bar;x2bar];
ubar=(-0.2*x1bar-0.3*x2bar)/0.2;

n = size(A,2); % tamanho de x
m = size(B,2); % tamanho de u
e = 0.001;

% P inicial
P = eye(n); 
i=0;
i_max=15;

while i<i_max
    % solving for P fixed
    prog = sosprogram(vars);
    [prog,K] = sospolymatrixvar(prog,monomials(x,0),[m,n]);
    % sos constraints
    M = A'*P + K'*B'*P + P*A + P*B*K;
    estabilidade = [xbar'*M*xbar xbar'*M; M*xbar M]+e*eye(n+1);
    prog = sosineq(prog, -estabilidade);
    % first solution
    sol = sossolve(prog);
    K = double(sosgetsol(sol,K));
    % solving for K fixed
    prog = sosprogram(vars);
    [prog,P] = sospolymatrixvar(prog,monomials(x,0),[n,n],'symmetric');
    % sos constraints
    prog = sosineq(prog, P-e*eye(n));
    M = A'*P + K'*B'*P + P*A + P*B*K;
    estabilidade = [xbar'*M*xbar xbar'*M; M*xbar M]+e*eye(n+1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
    prog = sosineq(prog, -estabilidade);
    % first solution
    sol = sossolve(prog);
    P = double(sosgetsol(sol,P));
    i=i+1;
end


pvar x1 x2
x = [x1;x2];
V = x'*P*x

%% malha aberta

syms x1 x2;
a = 4;
[x1,x2] = meshgrid(-a:0.01:a,-a:0.01:a);
dx1 = -0.5*(17.76 - 103.79.*x1 + 229.62.*x1.^2 - 226.31.*x1.^3 + 83.72.*x1.^4).*x1 + A(1,2).*x2;
dx2 = A(2,1).*x1 + A(2,2).*x2 + 0.24;
figure
streamslice(x1,x2,dx1,dx2,2)
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


% ts = 1e-11;
% tsim = 1e-9;
% t = 0:ts:tsim;
% x = zeros(n,length(t));
% x(:,1) = rand(2,1);
% for i=2:tsim
%     x(:,i) = x(:,i-1)+ts*[-0.5*(17.76 - 103.79*x(1,i-1) + 229.62*x(1,i-1)^2 - 226.31*x(1,i-1)^3 + 83.72*x(1,i-1)^4) 0.5;
%                           -0.2 -0.3]*x(:,i-1);
% end
% figure
% plot(t,x(1,:),t,x(2,:))
% legend('x1','x2','Location','best')




