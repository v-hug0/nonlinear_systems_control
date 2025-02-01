%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mobile robot control on a reference path 
%%% Victor Hugo, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Coeficientes do controlador 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pvar e1 e2 e4 e5
e_n=[e1;e2;e4;e5];
Ur = [0.1 0];
f = [Ur(1)*(e5+1); Ur(1)*e4; Ur(2)*(e5+1); -Ur(2)*e4];
g = [-1 e2; 0 -e1; 0 -e5-1; 0 e4];

n = size(f,2);
m = size(g,2);

% V inicial 

%z = monomials(e_n,1:2);
%V = z'*z;
V = (e_n'*e_n)^3
i = 0;
i_max = 100;
ubar = [0.1;0];
gamma = 10;
gamma_i = gamma;
while i<i_max
    % solving with SOSTOOLS (V,gamma fixed)
    prog = sosprogram(e_n);
    [prog,u1] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    [prog,u2] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    u=[u1;u2];
    [prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,l3] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)],'wscoeff');
    g1 = e4^2+e5^2-1;
    % SOS constraints
    gradV = [diff(V,e1) diff(V,e2) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*(u+ubar))-s1*g1-l2+gamma*l3;
    prog = sosineq(prog,dV);
    % first solution (u,l3)
    sol = sossolve(prog);
    u = sosgetsol(sol,u);
    l3 = sosgetsol(sol,l3);
    % solving with SOSTOOLS (u,l3 fixed)
    dpvar gamma;
    prog = sosprogram(e_n);
    prog = sosdecvar(prog,gamma);
    [prog,V] = sospolyvar(prog,[monomials(e_n,2:6)],'wscoeff');
    [prog,l1] = sospolyvar(prog,[monomials(e_n,2:6)],'wscoeff');
    [prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)], 'wscoeff');
    g1 = e4^2+e5^2-1;
    % SOS constraints
    prog = sosineq(prog,V);
    gradV = [diff(V,e1) diff(V,e2) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*(u+ubar))-s1*g1-l2+gamma*l3;
    prog = sosineq(prog,dV);
    prog = sosineq(prog,V-l1);
    prog = sosineq(prog,gamma_i-gamma);
    % first solution (V)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    gamma = sosgetsol(sol,gamma)
    gamma_i = gamma;
    if double(gamma)<0
        break
    end
    i = i+1
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
    
    
    v(k)=1.1912E-06*E1^4 + 2.1176E-08*E1^3*E2 - 5.1596E-10*E1^3*E4 + 2.7138E-06*E1^3*(E5-1) + 0.00010968*E1^2*E2^2 + 0.010223*E1^2*E2*E4 - 1.0829E-09*E1^2*E2*(E5-1) - 0.00025238*E1^2*E4^2 + 3.8988E-10*E1^2*E4*(E5-1) - 1.5641E-06*E1^2*(E5-1)^2 - 1.4439E-08*E1*E2^3 + 1.4534E-09*E1*E2^2*E4 - 1.6866E-05*E1*E2^2*(E5-1) - 3.9975E-09*E1*E2*E4^2 - 0.00059379*E1*E2*E4*(E5-1) + 9.2753E-10*E1*E2*(E5-1)^2 - 8.2093E-11*E1*E4^3 + 0.00058798*E1*E4^2*(E5-1) + 1.3136E-09*E1*E4*(E5-1)^2 + 2.4727E-06*E1*(E5-1)^3 - 7.4852E-05*E2^4 - 0.0069592*E2^3*E4 - 8.7049E-11*E2^3*(E5-1) + 0.00085724*E2^2*E4^2 + 2.5662E-09*E2^2*E4*(E5-1) + 5.0686E-06*E2^2*(E5-1)^2 - 0.00080217*E2*E4^3 - 4.6678E-10*E2*E4^2*(E5-1) + 0.00045399*E2*E4*(E5-1)^2 + 5.2068E-11*E2*(E5-1)^3 - 0.00027712*E4^4 - 4.0533E-09*E4^3*(E5-1) + 0.00086336*E4^2*(E5-1)^2 + 1.7393E-09*E4*(E5-1)^3 + 1.4806E-07*(E5-1)^4 + 0.097526*E1^3 + 5.0641E-08*E1^2*E2 - 1.4563E-07*E1^2*E4 + 0.0013808*E1^2*(E5-1) + 0.083149*E1*E2^2 - 0.022719*E1*E2*E4 + 6.3207E-08*E1*E2*(E5-1) + 0.030064*E1*E4^2 + 1.4549E-07*E1*E4*(E5-1) + 0.1448*E1*(E5-1)^2 - 3.1867E-09*E2^3 - 3.8842E-07*E2^2*E4 + 0.00049041*E2^2*(E5-1) - 7.9189E-09*E2*E4^2 - 0.0047597*E2*E4*(E5-1) - 2.9484E-09*E2*(E5-1)^2 - 4.7006E-07*E4^3 + 0.033397*E4^2*(E5-1) - 5.3557E-07*E4*(E5-1)^2 - 0.0012166*(E5-1)^3 + 0.00066429*E1^2 - 5.1593E-08*E1*E2 + 5.546E-08*E1*E4 + 0.019711*E1*(E5-1) - 0.0068108*E2^2 - 0.020545*E2*E4 - 3.4298E-08*E2*(E5-1) + 0.00036192*E4^2 - 2.2804E-08*E4*(E5-1) + 0.0088234*(E5-1)^2 + 0.092465*E1 - 1.0297E-07*E2 - 2.2051E-07*E4 + 0.099539*(E5-1) + ubar(1);
    w(k)=5.0288E-11*E1^4 - 0.011536*E1^3*E2 + 0.00028167*E1^3*E4 + 9.5628E-10*E1^3*(E5-1) + 3.4615E-10*E1^2*E2^2 + 1.4821E-09*E1^2*E2*E4 + 0.00075568*E1^2*E2*(E5-1) - 6.5765E-10*E1^2*E4^2 - 0.00068596*E1^2*E4*(E5-1) + 5.6368E-10*E1^2*(E5-1)^2 + 0.0078619*E1*E2^3 - 0.0009741*E1*E2^2*E4 - 2.974E-09*E1*E2^2*(E5-1) + 0.00090014*E1*E2*E4^2 + 5.0295E-10*E1*E2*E4*(E5-1) - 0.00050948*E1*E2*(E5-1)^2 + 0.00031282*E1*E4^3 + 4.5152E-09*E1*E4^2*(E5-1) - 0.00096954*E1*E4*(E5-1)^2 - 1.9221E-09*E1*(E5-1)^3 - 3.3827E-12*E2^4 + 7.7909E-10*E2^3*E4 - 4.3436E-05*E2^3*(E5-1) - 1.7196E-11*E2^2*E4^2 + 4.2637E-05*E2^2*E4*(E5-1) - 1.3364E-10*E2^2*(E5-1)^2 + 1.1813E-10*E2*E4^3 - 1.2931E-05*E2*E4^2*(E5-1) - 8.3771E-11*E2*E4*(E5-1)^2 - 6.5126E-06*E2*(E5-1)^3 + 6.4662E-11*E4^4 - 1.8719E-06*E4^3*(E5-1) - 1.2875E-10*E4^2*(E5-1)^2 + 4.3106E-06*E4*(E5-1)^3 + 4.8613E-11*(E5-1)^4 + 4.7018E-07*E1^3 + 0.041827*E1^2*E2 + 0.015515*E1^2*E4 - 3.0847E-07*E1^2*(E5-1) + 5.1413E-07*E1*E2^2 + 8.4259E-08*E1*E2*E4 + 0.0092799*E1*E2*(E5-1) + 1.8188E-07*E1*E4^2 - 0.025577*E1*E4*(E5-1) + 4.6421E-07*E1*(E5-1)^2 + 0.001827*E2^3 + 0.16709*E2^2*E4 + 1.2558E-09*E2^2*(E5-1) - 0.019286*E2*E4^2 - 1.424E-07*E2*E4*(E5-1) + 0.00078209*E2*(E5-1)^2 + 0.13451*E4^3 - 3.7289E-08*E4^2*(E5-1) + 0.1277*E4*(E5-1)^2 - 1.5505E-09*(E5-1)^3 - 8.0274E-08*E1^2 + 0.039743*E1*E2 + 0.0007619*E1*E4 + 1.6522E-07*E1*(E5-1) - 3.299E-08*E2^2 - 3.2336E-07*E2*E4 + 0.090801*E2*(E5-1) + 3.135E-09*E4^2 + 0.05408*E4*(E5-1) - 2.6656E-07*(E5-1)^2 + 2.1787E-07*E1 + 0.15487*E2 + 0.066543*E4 - 1.4074E-07*(E5-1) + ubar(2);
    

    X = [X x+Ts*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
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