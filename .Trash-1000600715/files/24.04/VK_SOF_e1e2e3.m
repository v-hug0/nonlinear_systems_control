clc
clear all

% This algorithm represents a SOF controller found by V-K iterative method, while optimizing a strictness multiplier.

pvar e1 e2 e3 e4 e5;
e=[e1;e2;e3;e4;e5];
Ur = [0.2 0.2];
 
f = [Ur(1)*(e5+1);Ur(1)*e4;Ur(2);Ur(2)*(e5+1);-Ur(2)*e4];
g = [-1 e2; 0 e1; 0 -1; 0 -e5-1; 0 e4];

% V inicial 

%V = beta*e'*e;
P = eye(length(e));
%V = e'*P*e;
V = 0.5*e'*P*e;
rho =  1;
ubar = [0.2;0.2];

nv = 1;
nt = 1;

%%
i = 0;
i_max = 5;
while i<i_max
    i
    %%%%%% solving for (V) fixed
    prog = sosprogram(e);
    %.....controller
    [prog,u1] = sospolyvar(prog,[monomials([e1;e2;e3],1:2)],'wscoeff');
    [prog,u2] = sospolyvar(prog,[monomials([e1;e2;e3],1:2)],'wscoeff');
    u=[u1;2];
    %.....S-procedure multipliers
    [prog,s1] = sospolyvar(prog, [monomials(e,1:4)],'wscoeff');
    [prog,alfa] = sospolyvar(prog,[monomials(e,2:6)],'wscoeff');
    %.....Strict term
    [prog,T] = sospolyvar(prog, [monomials(e,2:6)],'wscoeff');
    %.....SOS constraints
    g1 = e4^2+e5^2-1;  % estados amarrados
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    %...
    T_PD = T-1e-6*(e'*e)^nt;
    dV = -gradV*(f+g*(u+ubar))-alfa*(rho-V)-s1*g1-T; % dV after S-procedure
    %...
    prog = sosineq(prog,T_PD);
    prog = sosineq(prog,dV);
    % first solution (k,s1,alfa,T)
    sol = sossolve(prog);
    %...
    T = sosgetsol(sol,T);
    u = sosgetsol(sol,u);
    alfa = sosgetsol(sol,alfa);
    s1 = sosgetsol(sol,s1);
    %%%%%% solving for (k,T,alfa) fixed
    prog = sosprogram(e);
    %.....LF
    [prog,V] = sospolyvar(prog,[monomials(e,2:4)],'wscoeff');
    %.....S-procedure
    [prog,s1] = sospolyvar(prog, [monomials(e,1:4)],'wscoeff');
    %.....SOS constraints
    g1 = e4^2+e5^2-1;
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    %...
    V_PD = V-1e-6*(e'*e)^nv;
    dV = -gradV*(f+g*(u+ubar))-alfa*(rho-V)-s1*g1-T;
    %...
    prog = sosineq(prog,V_PD);
    prog = sosineq(prog,dV);
    % second solution (V,s1)
    sol = sossolve(prog);
    %...
    s1 = sosgetsol(sol,s1);
    V = sosgetsol(sol,V);
    i = i+1;
end





%%   

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Inicializacao de Variveis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 0.1;          % sampling time
Tsim = (2*pi*1)/Ts;   %90      % tempo de simula��o
h = 0.1;           % integrating time
D = 0.4;

%xr = [0.931386162550655; 1.388657897964; 1.98];  % x,y,theta referenciais para in�cio de trajet�ria. 
xr = [0;0;0];
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

xo = [0.2;0;pi/10];
x = xo;
X = x;
v = [];
w = [];
ubar = [0.2;0.2];

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

    v(k) = 0.0020794*E1^2 + 0.0024624*E1*E2 - 0.01153*E1*E3 - 0.0092007*E1*E4 - 0.0095686*E1*(E5-1) + 0.13263*E2^2 + 0.3332*E2*E3 + 0.72671*E2*E4 - 0.0020522*E2*(E5-1) + 0.011446*E3^2 + 0.0061094*E3*E4 - 3.0227E-05*E3*(E5-1) - 0.0058096*E4^2 - 2.1017E-05*E4*(E5-1) + 0.0011586*(E5-1)^2 + 0.74929*E1 + 0.38342*E2 + 0.0011503*E3 + 0.0034194*E4 + 0.1984*(E5-1) + ubar(1);
    w(k) = 0.024112*E1^2 + 0.022717*E1*E2 - 0.012547*E1*E3 + 0.014751*E1*E4 - 0.0001103*E1*(E5-1) - 0.024677*E2^2 - 0.00015734*E2*E3 - 0.0001721*E2*E4 + 0.00042034*E2*(E5-1) + 0.00073417*E3^2 - 0.0016208*E3*E4 - 0.0057224*E3*(E5-1) + 0.00020463*E4^2 - 0.043451*E4*(E5-1) + 0.0016931*(E5-1)^2 + 0.00041107*E1 + 0.07407*E2 + 0.19685*E3 + 0.42262*E4 - 8.1994E-05*(E5-1) + ubar(2);


    X = [X x+r(k)+Ts*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
    
    Erro(:,k) = Xr(:,k+1)- X(:,k+1);
    Erro4(k) = E4;
    Erro5(k) = E5;
   
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
subplot(2,1,1)
plot(tempo,Erro(1,:),tempo,Erro(2,:),tempo,Erro(3,:));
legend('e1','e2','e3')

hold on
subplot(2,1,2)
plot(tempo,Erro4,tempo,Erro5)
legend('e4','e5')

 EE1 = Erro(1,:);
 EE2 = Erro(2,:);
 EE3 = Erro(3,:);





