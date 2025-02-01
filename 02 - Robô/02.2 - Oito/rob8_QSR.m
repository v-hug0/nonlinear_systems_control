%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mobile robot control on a reference path 
%%% Gregor Klancar, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modificado por Marcus
% 08.04.2016

clear all
close all
clc
%%
format long g     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Inicializacao de Variveis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 0.1;          % sampling time
Tsim = (2*pi*1)/Ts;   %90      % tempo de simula��o
h = 0.1;           % integrating time
D = 0.4;

xr = [0; 0; 0];  % x,y,theta referenciais para in�cio de trajet�ria. 
Xr = xr;         %Euler

Tsim = 85;
% The next loop difines the reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Geracao da trajetoria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:round(Tsim/h)  % 1 a tempo de simula��o/ tempo de integra��o   
    vr(k) = 0.3;
    if k<= 419
        wr(k) = 0.15;
     else
         wr(k) = -0.15;
    end
    Xr = [Xr xr+h*([vr(k)*cos(xr(3));vr(k)*sin(xr(3));wr(k)])];  %Euler
    xr = Xr(:,k+1);
end
Ur = [vr;wr];
hold on
plot(Xr(1,:),Xr(2,:),'')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Coeficientes do controlador 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pvar e1 e2 e3 u1 u2
var = [e1;e2;e3;u1;u2];
e = [e1;e2;e3];
u = [u1;u2];
Ur = [0.2  0.2];

f = [Ur(1)*(1-0.5*e3^2+0.042*e3^4-0.0014*e3^6);
     Ur(1)*(e3-0.17*e3^3+0.0083*e3^5-0.0002*e3^7);
     Ur(2)];
g = [-1 e2; 0 -e1; 0 -1];

h = monomials(e,1:2);

n = size(f,2);
m = size(g,2);
p = size(h,1); %numero de linhas do h
j = size(h,2);

betaV = 1e-4;
betaT = 1e-4;
betaL = 1e-6;
nv = 2;
nt = 2;
nl = 2;
alfa0 = 1e-6*(e1^2+e2^2+e3^2);

i = 0;
imax = 100; 

prog=sosprogram(var);
[prog,V0] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
[prog,T0] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
%[prog,L0] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
[prog,Q0]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
[prog,S0]=sospolymatrixvar(prog,monomials(var,0),[p,m]);
[prog,R0]=sospolymatrixvar(prog,monomials(var,0),[m,m],'symmetric');

gradV0 = [diff(V0,e1) diff(V0,e2) diff(V0,e3)];
%cond1 = -gradV0*(f+g*u)-T0+h'*Q0*h+2*h'*S0*u+u'*R0*u-alfa0*(1-L0);
cond1 = -gradV0*(f+g*u)-T0+h'*Q0*h+2*h'*S0*u+u'*R0*u;
cond2 = V0-betaV*(e1^2+e2^2+e3^2)^nv;
cond3 = T0-betaT*(e1^2+e2^2+e3^2)^nt;
%cond4 = L0-betaL*(e1^2+e2^2+e3^2)^nl;
prog = sosineq(prog, cond1);
prog = sosineq(prog, cond2);
prog = sosineq(prog, cond3);
%prog = sosineq(prog, cond4);
prog = sosineq(prog, R0);

sol = sossolve(prog);

Q0 = double(sosgetsol(sol,Q0));
R0 = double(sosgetsol(sol,R0));
S0 = double(sosgetsol(sol,S0));
%L0 = sosgetsol(sol,L0);

if S0*inv(R0)*S0'-Q0 >=0
    K = -inv(R0)*S0';
    V0 = sosgetsol(sol,V0);
else 
    while i<imax
    prog=sosprogram(var);
    [prog,V] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
    [prog,T] = sospolyvar(prog,monomials(e,2:8),'wscoeff');
    %[prog,alfa] = sospolyvar(prog,monomials(e,2:8),'wscoeff');
    [prog,Q]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
    [prog,S]=sospolymatrixvar(prog,monomials(var,0),[p,m]);
    [prog,R]=sospolymatrixvar(prog,monomials(var,0),[m,m],'symmetric');
    
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3)];
    %cond1 = -gradV*(f+g*u)-T+h'*Q*h+2*h'*S*u+u'*R*u-alfa*(1-L0);
    cond1 = -gradV*(f+g*u)-T+h'*Q*h+2*h'*S*u+u'*R*u;
    cond2 = V-betaV*(e1^2+e2^2+e3^2)^nv;
    cond3 = T-betaT*(e1^2+e2^2+e3^2)^nt;
    cond5 = R0-R;
    cond6 = S*inv(R0)*S0'+S0*inv(R0)*S'-2*S0*inv(R0)*S0'+Q0-Q;
    prog = sosineq(prog, cond1);
    prog = sosineq(prog, cond2);
    prog = sosineq(prog, R);
    prog = sosineq(prog, cond5);
    prog = sosineq(prog, cond6);


    sol = sossolve(prog);

    Q = double(sosgetsol(sol,Q));
    R = double(sosgetsol(sol,R));
    S = double(sosgetsol(sol,S));
    %alfa = sosgetsol(sol,alfa);

    Q0 = Q;
    S0 = S;
    R0 = R;
    %alfa0 = alfa;

    d = S*inv(R)*S'-Q;
    min(eig(d))
    i=i+1

    if min(eig(d)) >= 0 
        K = -inv(R)*S';
        V = sosgetsol(sol,V);
        break
    end
end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Simulacao do robo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [-0.5;0.5;0];
x = xo;
X = x;
v = [];
w = [];
ubar = [0.3;0.15];
ubar2 = [0.3;-0.15];

m = 0;
d = 0.000;
r = m + d.*randn(1,round(Tsim/h));

k1 = 2.6;
k2 = 1;
k3 = 3;

for k = 1:round(Tsim/h)    % 1 a tempo de simula��o/ tempo de integra��o  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(1);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);
    
    
    if k<= 500
        v(k) = -0.020288*E1^2 + 0.010065*E1*E2 + 0.019745*E1*E3 - 0.01399*E2^2 + 0.017534*E2*E3 + 0.0068482*E3^2 + 0.23245*E1 + 0.011004*E2 - 0.26535*E3 + ubar(1);
        w(k) = 0.01516*E1^2 + 0.29945*E1*E2 + 0.028598*E1*E3 - 0.011502*E2^2 + 0.20065*E2*E3 - 0.064419*E3^2 - 0.23641*E1 - 0.010334*E2 + 0.57196*E3 + ubar(2);
    else
        %v(k) =  + ubar2(1);
        %w(k) =  + ubar2(2);
        %v(k) = -0.019638*E1^2 + 0.0025691*E1*E2 + 0.019375*E1*E3 - 0.016991*E2^2 + 0.0032373*E2*E3 + 0.00076959*E3^2 + 0.17215*E1 + 0.0093107*E2 - 0.1612*E3 + ubar2(2);
        %w(k) = 0.014328*E1^2 + 0.20234*E1*E2 + 0.005553*E1*E3 - 0.007189*E2^2 + 0.11527*E2*E3 - 0.044053*E3^2 - 0.12835*E1 - 0.006847*E2 + 0.32723*E3 + ubar2(2);

    end

    if k>550 && k<555
        X = [X x+r(k)+Ts*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    else
        X = [X x+r(k)+Ts*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    end

    x = X(:,k+1);
    
    Erro(:,k) = Xr(:,k+1)- X(:,k+1);
   
end


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
