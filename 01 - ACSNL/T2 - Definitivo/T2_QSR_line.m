clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Coeficientes do controlador 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pvar e1 e2 e4 e5 u1 u2
var = [e1;e2;e4;e5;u1;u2];
e = [e1;e2;e4;e5];
u = [u1;u2];
Ur = [0.2 0.2];

f = [Ur(1)*(e5+1);
     Ur(1)*e4;
     Ur(2)*(e5+1);
    -Ur(2)*e4];
g = [-1 e2;
      0 e1;
      0 -e5-1;
      0 e4];

%h = monomials(e,1:3);
h1 = monomials(e,1);
h1 = ones(1,length(h1))*h1;
h2 = monomials(e,2);
h1 = ones(1,length(h2))*h2;
h3 = monomials(e,3);
h3 = ones(1,length(h3))*h3;
h4 = monomials(e,4);
h4 = ones(1,length(h4))*h4;
h = [h1;h2;h3;h4];

n = size(f,2);
m = size(g,2);
p = size(h,1);

ubar = [0.2;0.2];
i = 0;
imax = 100; 

S0 = zeros(p,m);
R0 = eye(m);
Ls = [-inv(R0)*S0' -eye(m,m)]';
lambda = 10;
lambda_i = lambda;

while i<imax
    dpvar lambda
    prog=sosprogram(var);
    prog = sosdecvar(prog,lambda);
    [prog,V] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
    [prog,l1] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
    [prog,l2] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
    [prog,s1] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
    
    [prog,Q]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
    [prog,S]=sospolymatrixvar(prog,monomials(var,0),[p,m]);
    [prog,R]=sospolymatrixvar(prog,monomials(var,0),[m,m]);
    %R = [R1(1) 0; 0 R1(2)];
    
    g1 = e4^2+e5^2-1;
    gradV = [diff(V,e1) diff(V,e2) diff(V,e4) diff(V,e5)];
    cond1 = V-l1;
    cond2 = gradV*(f+g*(u+ubar))-h'*Q*h-2*h'*S*(u+ubar)-(u+ubar)'*R*(u+ubar)+s1*g1+l2;
    cond3 = [Q S ;S' R] + Ls*[S' R] + [S;R]*Ls'+lambda*[-eye(p) zeros(p,m); zeros(m,p) zeros(m,m)];
    prog = sosineq(prog, cond1);
    prog = sosineq(prog, -cond2);
    prog = sosineq(prog, -cond3);
    prog = sosineq(prog, lambda_i-lambda);
    
    sol = sossolve(prog);
    
    Q = double(sosgetsol(sol,Q));
    R = double(sosgetsol(sol,R));
    S = double(sosgetsol(sol,S));

    lambda = double(sosgetsol(sol,lambda))
    lambda_i = lambda;
    
    Q0 = Q;
    S0 = S;
    R0 = R;
    Ls = [-inv(R0)*S0' -eye(m,m)]';
    
    d = S*inv(R)*S'-Q;
    min(eig(d))
    i=i+1
    
    if lambda<0 || min(eig(d)) >= 0 
    %if min(eig(d)) >= 0 
        K = -inv(R)*S';
        V = sosgetsol(sol,V);
        break
    end
end


%%   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Inicializacao de Variveis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 0.1;          % sampling time
Tsim = (2*pi*1)/Ts;   %90      % tempo de simula��o
h = 0.1;           % integrating time
D = 0.4;

xr = [2.5; 0.3; 0];  % x,y,theta referenciais para in�cio de trajet�ria. 
Xr = xr;         %Euler

Tsim = 200;
lado = 2;
% The next loop difines the reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Geracao da trajetoria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:round(Tsim/h)  % 1 a tempo de simula��o/ tempo de integra��o   
    vr(k)=0.1;
    wr(k)=0;

    Xr = [Xr xr+h*([vr(k)*cos(xr(3));vr(k)*sin(xr(3));wr(k)])];  %Euler
    xr = Xr(:,k+1);
end
Ur = [vr;wr];
hold on
plot(Xr(1,:),Xr(2,:),'')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Simulacao do robo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [0;0;0];
x = xo;
X = x;
v = [];
w = [];
ubar = [0.1;0];
vr = 0.1;

for k = 1:round(Tsim/h)    % 1 a tempo de simula��o/ tempo de integra��o  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(1);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);
    
    
    v(k) = 0.58036*E1^3 - 7.5167E-06*E1^2*E2 + 9.833E-08*E1^2*E4 + 0.0021413*E1^2*(E5-1) - 0.3743*E1*E2^2 - 0.012165*E1*E2*E4 + 5.3854E-07*E1*E2*(E5-1) + 0.052015*E1*E4^2 + 3.1334E-07*E1*E4*(E5-1) + 0.058978*E1*(E5-1)^2 - 2.6572E-06*E2^3 + 9.5191E-07*E2^2*E4 + 0.0038964*E2^2*(E5-1) - 3.735E-07*E2*E4^2 - 0.00039372*E2*E4*(E5-1) - 7.4755E-07*E2*(E5-1)^2 - 2.9763E-06*E4^3 + 0.00029775*E4^2*(E5-1) - 8.7714E-07*E4*(E5-1)^2 + 0.0013213*(E5-1)^3 - 0.06266*E1^2 - 9.913E-07*E1*E2 + 2.3025E-06*E1*E4 + 0.022247*E1*(E5-1) - 0.3955*E2^2 + 0.022489*E2*E4 - 5.2033E-07*E2*(E5-1) - 0.13653*E4^2 + 7.0985E-07*E4*(E5-1) - 0.17728*(E5-1)^2 + 0.046396*E1 - 5.591E-07*E2 - 2.9131E-06*E4 + 0.0086738*(E5-1) + ubar(1);
    w(k) = -2.8665E-06*E1^3 + 1.0583*E1^2*E2 + 0.077856*E1^2*E4 - 8.905E-07*E1^2*(E5-1) + 1.4536E-06*E1*E2^2 - 5.2826E-07*E1*E2*E4 + 0.01688*E1*E2*(E5-1) - 3.8989E-07*E1*E4^2 + 0.070336*E1*E4*(E5-1) + 9.4341E-07*E1*(E5-1)^2 + 0.68256*E2^3 + 0.0033938*E2^2*E4 - 9.6394E-07*E2^2*(E5-1) + 0.17707*E2*E4^2 + 1.0455E-06*E2*E4*(E5-1) + 0.23842*E2*(E5-1)^2 + 0.37132*E4^3 - 2.1059E-07*E4^2*(E5-1) + 0.17643*E4*(E5-1)^2 + 5.1424E-07*(E5-1)^3 - 6.6129E-07*E1^2 - 0.012385*E1*E2 - 0.1943*E1*E4 + 8.9214E-08*E1*(E5-1) + 2.0154E-06*E2^2 - 1.8953E-07*E2*E4 + 0.038586*E2*(E5-1) + 1.5839E-06*E4^2 - 0.009757*E4*(E5-1) + 6.0472E-07*(E5-1)^2 - 1.4609E-06*E1 + 0.19104*E2 + 0.33663*E4 - 1.3074E-06*(E5-1) + ubar(2);
    

    X = [X x+Ts*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
    
    Erro(:,k) = Xr(:,k+1)- X(:,k+1);
   
end


figure
plot(Xr(1,:),Xr(2,:),'')
hold on
plot(X(1,:),X(2,:),'r')
legend('rungeKutta','euler','trajet�ria do rob�')
xlabel('x[m]')
ylabel('y[m]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem de v e w 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
ne = length(Erro);
tempo = Ts*(0:ne-1);
subplot(2,1,1)
plot(tempo,v,'b')
legend('v')

hold on
subplot(2,1,2)
plot (tempo,w,'r')
legend('w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem do erro x y e theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(tempo,Erro(1,:),tempo,Erro(2,:),tempo,Erro(3,:));
legend('x','y','theta')

 EE1 = Erro(1,:);
 EE2 = Erro(2,:);
 EE3 = Erro(3,:);