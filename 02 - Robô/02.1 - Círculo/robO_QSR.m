%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mobile robot control on a reference path 
%%% Gregor Klancar, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modificado por Marcus
% 08.04.2016

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

h1 = monomials(e,1);
h1 = ones(1,length(h1))*h1;
h2 = monomials(e,2);
h2 = ones(1,length(h2))*h2;
h3 = monomials(e,3);
h3 = ones(1,length(h3))*h3;
h4 = monomials(e,4);
h4 = ones(1,length(h4))*h4;

h = [h1;h2;h3;h4];

n = size(f,2);
m = size(g,2);
p = size(h,1); %numero de linhas do h
j = size(h,2);

ubar = [0.2;0.2];
i = 0;
imax = 100; 

prog=sosprogram(var);
[prog,V0] = sospolyvar(prog,monomials(e,2),'wscoeff');
[prog,l1] = sospolyvar(prog,monomials(e,2),'wscoeff');
[prog,l2] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
[prog,s1] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
%[prog,L0] = sospolyvar(prog,monomials(e,2:4),'wscoeff');

[prog,Q0]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
[prog,S]=sospolymatrixvar(prog,monomials(var,0),[p,m]);
% [prog,S1]=sospolymatrixvar(prog,monomials(var,0),[20,1]); 
% [prog,S2]=sospolymatrixvar(prog,monomials(var,0),[23,1]);
% S = [[S1(1);S1(2);S1(3);S1(4); 
%       S1(5);S1(6);S1(7);0;0;0;0;0;0;S1(8);
%       S1(9);S1(10);S1(11);0;S1(12);S1(13);S1(14);S1(15);S1(16);0;0;0;S1(17);0;S1(18);S1(19);S1(20);0;0;0] [S2(1);S2(2);S2(3);S2(4);  
%       0;S2(5);0;S2(6);0;0;0;0;0;S2(7);
%       S2(8);S2(9);0;S2(10);S2(11);0;S2(12);S2(13);S2(14);0;S2(15);S2(16);S2(17);S2(18);S2(19);S2(20);0;S2(21);S2(22);S2(23)]]; 
[prog,R1]=sospolymatrixvar(prog,monomials(var,0),[m,m]); 
R = [R1(1) 0; 0 R1(2)]; 



g1 = e4^2+e5^2-1;
gradV0 = [diff(V0,e1) diff(V0,e2) diff(V0,e4) diff(V0,e5)];
cond1 = -gradV0*(f+g*(u))+h'*Q0*h+2*h'*S*u+u'*R*u-s1*g1-l2;
%cond1 = -gradV0*(f+g*(u+ubar))+h'*Q0*h+2*h'*S0*(u+ubar)+(u+ubar)'*R0*(u+ubar)-alfa0*(1-L0)-s1*g1-l2;
cond2 = V0-l1;
%cond4 = L0-betaL*(e'*e)^nl;
prog = sosineq(prog, cond1);
prog = sosineq(prog, cond2);
%prog = sosineq(prog, cond3);
%prog = sosineq(prog, cond4);
prog = sosineq(prog, R);

sol = sossolve(prog);

Q0 = double(sosgetsol(sol,Q0));
R0 = double(sosgetsol(sol,R));
S0 = double(sosgetsol(sol,S));
%L0 = sosgetsol(sol,L0);


if S0*inv(R0)*S0'-Q0 >=0
    K = -inv(R0)*S0';
    V0 = sosgetsol(sol,V0);
else 
    while i<imax
    prog=sosprogram(var);
    [prog,V] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
    [prog,l1] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
    [prog,l2] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
    [prog,s1] = sospolyvar(prog,monomials(e,2:4),'wscoeff');
    %[prog,alfa] = sospolyvar(prog,monomials(e,2:6),'wscoeff');

    [prog,Q]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
    [prog,S]=sospolymatrixvar(prog,monomials(var,0),[p,m]);
    % [prog,S1]=sospolymatrixvar(prog,monomials(var,0),[20,1]); 
    % [prog,S2]=sospolymatrixvar(prog,monomials(var,0),[23,1]);
    % S = [[S1(1);S1(2);S1(3);S1(4); 
    %       S1(5);S1(6);S1(7);0;0;0;0;0;0;S1(8);
    %       S1(9);S1(10);S1(11);0;S1(12);S1(13);S1(14);S1(15);S1(16);0;0;0;S1(17);0;S1(18);S1(19);S1(20);0;0;0] [S2(1);S2(2);S2(3);S2(4);  
    %       0;S2(5);0;S2(6);0;0;0;0;0;S2(7);
    %       S2(8);S2(9);0;S2(10);S2(11);0;S2(12);S2(13);S2(14);0;S2(15);S2(16);S2(17);S2(18);S2(19);S2(20);0;S2(21);S2(22);S2(23)]]; 
    % [prog,R1]=sospolymatrixvar(prog,monomials(var,0),[m,m]); 
    [prog,R]=sospolymatrixvar(prog,monomials(var,0),[m,m],'symmetric');

    g1 = e4^2+e5^2-1;
    gradV = [diff(V,e1) diff(V,e2) diff(V,e4) diff(V,e5)];
    %cond1 = -gradV*(f+g*(u+ubar))+h'*Q*h+2*h'*S*(u+ubar)+(u+ubar)'*R*(u+ubar)-alfa*(1-L0)-s1*g1-l2;
    cond1 = -gradV*(f+g*(u))+h'*Q*h+2*h'*S*(u)+(u)'*R*(u)-s1*g1-l2;
    cond2 = V-l1;
    %cond2 = V-betaV*(e'*e)^nv;
    %cond3 = l2-betaT*(e'*e)^nt;
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
    d_i = d;
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
format long g     

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Inicializacao de Variveis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 0.1;          % sampling time
Tsim = (2*pi*1)/Ts;   %90      % tempo de simula��o
h = 0.1;           % integrating time
D = 0.4;

%xr = [0.931386162550655; 1.388657897964; 1.98];  % x,y,theta referenciais para in�cio de trajet�ria. 
xr = [1; 1; 0];
Xr = xr;         %Euler

% The next loop difines the reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Geracao da trajetoria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:round(Tsim/h)  % 1 a tempo de simula��o/ tempo de integra��o   
    vr(k) = 0.2;
    wr(k) = 0.2;
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
ubar = [0.2;0.2];

%limitV = 0.3;
%limitW = 0.3;

%[MIN,ki]=min(abs(Xr(1,:)-xo(1)));


m = 0;
d = 0.000;
r = m + d.*randn(1,round(Tsim/h));

for k = 1:(round(Tsim/h))    % 1 a tempo de simula��o/ tempo de integra��o  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(1);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);
    
    v(k)=-0.027808*E1^2 - 0.075837*E1*E2 + 0.1456*E1*E4 - 0.074339*E1*(E5-1) + 0.0064466*E2^2 + 0.077405*E2*E4 - 0.044482*E2*(E5-1) - 0.20916*E4^2 - 0.052158*E4*(E5-1) - 0.13144*(E5-1)^2 + 1.2835*E1 + 0.55863*E2 - 1.3302*E4 - 0.43386*(E5-1) + ubar(1);
    w(k)=-0.43998*E1^2 - 2.3658*E1*E2 + 0.44347*E1*E4 - 1.1488*E1*(E5-1) - 0.58879*E2^2 + 0.85038*E2*E4 - 0.1447*E2*(E5-1) - 0.79959*E4^2 - 1.1407*E4*(E5-1) + 0.027343*(E5-1)^2 - 0.971*E1 - 0.36678*E2 + 1.0347*E4 + 0.33595*(E5-1) + ubar(2);


    v(k) = 1.4629*E1^4 + 1.4629*E1^3*E2 + 1.4629*E1^3*E4 + 1.4629*E1^3*(E5-1) + 1.4629*E1^2*E2^2 + 1.4629*E1^2*E2*E4 + 1.4629*E1^2*E2*(E5-1) + 1.4629*E1^2*E4^2 + 1.4629*E1^2*E4*(E5-1) + 1.4629*E1^2*(E5-1)^2 + 1.4629*E1*E2^3 + 1.4629*E1*E2^2*E4 + 1.4629*E1*E2^2*(E5-1) + 1.4629*E1*E2*E4^2 + 1.4629*E1*E2*E4*(E5-1) + 1.4629*E1*E2*(E5-1)^2 + 1.4629*E1*E4^3 + 1.4629*E1*E4^2*(E5-1) + 1.4629*E1*E4*(E5-1)^2 + 1.4629*E1*(E5-1)^3 + 1.4629*E2^4 + 1.4629*E2^3*E4 + 1.4629*E2^3*(E5-1) + 1.4629*E2^2*E4^2 + 1.4629*E2^2*E4*(E5-1) + 1.4629*E2^2*(E5-1)^2 + 1.4629*E2*E4^3 + 1.4629*E2*E4^2*(E5-1) + 1.4629*E2*E4*(E5-1)^2 + 1.4629*E2*(E5-1)^3 + 1.4629*E4^4 + 1.4629*E4^3*(E5-1) + 1.4629*E4^2*(E5-1)^2 + 1.4629*E4*(E5-1)^3 + 1.4629*(E5-1)^4 - 1.2459*E1^3 - 1.2459*E1^2*E2 - 1.2459*E1^2*E4 - 1.2459*E1^2*(E5-1) - 1.2459*E1*E2^2 - 1.2459*E1*E2*E4 - 1.2459*E1*E2*(E5-1) - 1.2459*E1*E4^2 - 1.2459*E1*E4*(E5-1) - 1.2459*E1*(E5-1)^2 - 1.2459*E2^3 - 1.2459*E2^2*E4 - 1.2459*E2^2*(E5-1) - 1.2459*E2*E4^2 - 1.2459*E2*E4*(E5-1) - 1.2459*E2*(E5-1)^2 - 1.2459*E4^3 - 1.2459*E4^2*(E5-1) - 1.2459*E4*(E5-1)^2 - 1.2459*(E5-1)^3 - 2.8433*E1^2 + 10.9321*E1*E2 - 8.4985*E1*E4 - 6.5765*E1*(E5-1) - 2.5089*E2^2 - 0.98483*E2*E4 + 2.1264*E2*(E5-1) - 3.233*E4^2 - 2.2591*E4*(E5-1) - 17.1039*(E5-1)^2 + ubar(1);
    w(k) = -1.084*E1^4 - 1.084*E1^3*E2 - 1.084*E1^3*E4 - 1.084*E1^3*(E5-1) - 1.084*E1^2*E2^2 - 1.084*E1^2*E2*E4 - 1.084*E1^2*E2*(E5-1) - 1.084*E1^2*E4^2 - 1.084*E1^2*E4*(E5-1) - 1.084*E1^2*(E5-1)^2 - 1.084*E1*E2^3 - 1.084*E1*E2^2*E4 - 1.084*E1*E2^2*(E5-1) - 1.084*E1*E2*E4^2 - 1.084*E1*E2*E4*(E5-1) - 1.084*E1*E2*(E5-1)^2 - 1.084*E1*E4^3 - 1.084*E1*E4^2*(E5-1) - 1.084*E1*E4*(E5-1)^2 - 1.084*E1*(E5-1)^3 - 1.084*E2^4 - 1.084*E2^3*E4 - 1.084*E2^3*(E5-1) - 1.084*E2^2*E4^2 - 1.084*E2^2*E4*(E5-1) - 1.084*E2^2*(E5-1)^2 - 1.084*E2*E4^3 - 1.084*E2*E4^2*(E5-1) - 1.084*E2*E4*(E5-1)^2 - 1.084*E2*(E5-1)^3 - 1.084*E4^4 - 1.084*E4^3*(E5-1) - 1.084*E4^2*(E5-1)^2 - 1.084*E4*(E5-1)^3 - 1.084*(E5-1)^4 - 0.98912*E1^3 - 0.98912*E1^2*E2 - 0.98912*E1^2*E4 - 0.98912*E1^2*(E5-1) - 0.98912*E1*E2^2 - 0.98912*E1*E2*E4 - 0.98912*E1*E2*(E5-1) - 0.98912*E1*E4^2 - 0.98912*E1*E4*(E5-1) - 0.98912*E1*(E5-1)^2 - 0.98912*E2^3 - 0.98912*E2^2*E4 - 0.98912*E2^2*(E5-1) - 0.98912*E2*E4^2 - 0.98912*E2*E4*(E5-1) - 0.98912*E2*(E5-1)^2 - 0.98912*E4^3 - 0.98912*E4^2*(E5-1) - 0.98912*E4*(E5-1)^2 - 0.98912*(E5-1)^3 - 4.5375*E1^2 - 30.3061*E1*E2 - 7.3743*E1*E4 - 4.0263*E1*(E5-1) - 5.4038*E2^2 - 3.9423*E2*E4 - 8.9445*E2*(E5-1) - 10.5131*E4^2 - 12.1807*E4*(E5-1) - 2.0376*(E5-1)^2 + ubar(2);

    % dissipatividade
    %v(k) = 0.031936*E1^2 - 0.028295*E1*E2 - 0.14718*E1*E3 - 0.014605*E2^2 - 0.066124*E2*E3 + 0.13842*E3^2 + 0.52793*E1 + 0.016323*E2 - 0.36578*E3 + ubar(1);
    %w(k) = -0.035067*E1^2 + 0.040549*E1*E2 + 0.232*E1*E3 - 0.040055*E2^2 + 0.24684*E2*E3 - 0.23467*E3^2 - 0.23177*E1 + 0.031981*E2 + 0.42755*E3 + ubar(2);
        
    % if v(k) >= limitV
    %     v(k) = limitV;
    % end
    % if w(k) >= limitW
    %     w(k) = limitW;
    % end
    % if v(k) < -limitV
    %     v(k) = -limitV;
    % end
    % if w(k) < -limitW
    %     w(k) = -limitW;
    % end


    X = [X x+r(k)+Ts*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
    
    Erro(:,k) = Xr(:,k+1)- X(:,k+1);
   
end


plot(X(1,:),X(2,:),'r')
legend('rungeKutta','trajetoria do robo')
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
