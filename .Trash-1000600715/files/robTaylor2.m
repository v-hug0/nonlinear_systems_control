%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mobile robot control on a reference path 
%%% Victor Hugo, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
%%
format long g     

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Inicializacao de Variveis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 0.1;          % sampling time
Tsim = (2*pi*1)/Ts;   %90      % tempo de simula��o
h = 0.1;           % integrating time
D = 0.4;

xr = [0; 0; 0];  % x,y,theta referenciais para in�cio de trajet�ria. 
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Coeficientes do controlador 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pvar e1 e2 e3 e4 e5 vr wr
e_n=[e1;e2;e3;e4;e5];
Ur = [0.3 -0.15];
%f = [Ur(1)*(e5+1); Ur(1)*e4; Ur(2); Ur(2)*(e5+1); -Ur(2)*e4];
f = [vr*(e5+1); vr*e4; wr; wr*(e5+1); -wr*e4];
g = [-1 e2; 0 -e1; 0 -1; 0 -e5-1; 0 e4];

n = size(f,2);
m = size(g,2);

% V inicial 
V = e_n'*e_n;
i = 0;
i_max = 20;
while i<i_max
    % solving with SOSTOOLS (V fixed)
    prog = sosprogram([e_n;vr;wr]);
    [prog,kx1] = sospolyvar(prog,[monomials(e_n,1:2)],'wscoeff');
    [prog,kx2] = sospolyvar(prog,[monomials([e2;e3;e4],1:2)],'wscoeff');
    kx=[kx1;kx2];
    [prog,T] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,1:4)],'wscoeff');
    g1 = e4^2+e5^2-1;
    % SOS constraints
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*kx)-T-s1*g1;
    prog = sosineq(prog,dV);
    prog = sosineq(prog,-vr+1);
    prog = sosineq(prog,vr+1);
    prog = sosineq(prog,-wr+1);
    prog = sosineq(prog,wr+1);
    % first solution (k and lambda)
    sol = sossolve(prog);
    kx = sosgetsol(sol,kx);
    % solving with SOSTOOLS (k fixed)
    prog = sosprogram(e_n);
    [prog,V] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,1:4)], 'wscoeff');
    [prog,T] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    g1 = e4^2+e5^2-1;
    % SOS constraints
    prog = sosineq(prog,V);
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*kx)-T-s1*g1;
    prog = sosineq(prog,dV);
    prog = sosineq(prog,-vr+1);
    prog = sosineq(prog,vr+1);
    prog = sosineq(prog,-wr+1);
    prog = sosineq(prog,wr+1);
    % first solution (V)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    i = i+1
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Simulacao do robo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [0.5;0.3;-0.4];
x = xo;
X = x;
v = [];
w = [];
ubar = [0.2;0.2];



%[MIN,ki]=min(abs(Xr(1,:)-xo(1)));


m = 0;
d = 0.000;
r = m + d.*randn(1,round(Tsim/h));

vr = 0.2;
wr = 0.2;
k1 = 2.6;
k2 = 1;
k3 = 3;

for k = 1:(round(Tsim/h))    % 1 a tempo de simula��o/ tempo de integra��o  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(1);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
   
    w(k) = wr+0.5*vr*(k2*(E2+k3*E3)+(1/k3)*(E3-0.17*E3^3+0.0083*E3^5-0.0002*E3^7));
    v(k) = vr*(1-0.5*E3^2+0.042*E3^4-0.0014*E3^6)-k3*E3*w(k)+k1*E1;

    %w(k) = 0.2+0.5*0.2*(3*E2+2*3*E3+0.33*(E3-0.17*E3^3+0.0083*E3^5-0.0002*E3^7));
    %v(k) = 0.2*(1-0.5*E3^2+0.042*E3^4-0.0014*E3^6)-E3*w(k)+E1;


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