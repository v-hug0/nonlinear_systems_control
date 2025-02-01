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
V = e_n'*e_n
i = 0;
i_max = 10;
ubar = [0.1;0];
gamma = 50;
while i<i_max
    % solving with SOSTOOLS (V fixed)
    %dpvar gamma;
    prog = sosprogram(e_n);
    %prog = sosdecvar(prog,gamma);
    [prog,u1] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    [prog,u2] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    u=[u1;u2];
    [prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)],'wscoeff');
    %[prog,s2] = sospolyvar(prog, [monomials(e_n,2:6)],'wscoeff');
    g1 = e4^2+e5^2-1;
    % SOS constraints
    gradV = [diff(V,e1) diff(V,e2) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*(u+ubar))-s1*g1-l2;
    prog = sosineq(prog,dV);
    %prog = sosineq(prog,gamma);
    %prog = sossetobj(prog,gamma);
    % first solution (k and s2)
    sol = sossolve(prog);
    u = sosgetsol(sol,u);
    %s2 = sosgetsol(sol,s2);
    %gamma = sosgetsol(sol,gamma)
    % solving with SOSTOOLS (k,s2 fixed)
    prog = sosprogram(e_n);
    [prog,V] = sospolyvar(prog,[monomials(e_n,2:6)],'wscoeff');
    [prog,l1] = sospolyvar(prog,[monomials(e_n,2:6)],'wscoeff');
    [prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)], 'wscoeff');
  
    g1 = e4^2+e5^2-1;
    % SOS constraints
    prog = sosineq(prog,V);
    gradV = [diff(V,e1) diff(V,e2) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*(u+ubar))-s1*g1-l2;
    prog = sosineq(prog,dV);
    prog = sosineq(prog,V-l1);
    %prog = sosineq(prog,T-l2);
    % first solution (V)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
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
    
    
    %v(k)= 0.0019098*E1^4 + 0.0015746*E1^3*E2 + 0.0073297*E1^3*E4 + 0.0020766*E1^3*(E5-1) - 0.005046*E1^2*E2^2 - 0.014478*E1^2*E2*E4 + 0.004803*E1^2*E2*(E5-1) - 0.0087539*E1^2*E4^2 + 0.011048*E1^2*E4*(E5-1) - 0.00042612*E1^2*(E5-1)^2 - 0.001336*E1*E2^3 - 0.0043726*E1*E2^2*E4 - 0.0033933*E1*E2^2*(E5-1) - 0.0015347*E1*E2*E4^2 - 0.0059516*E1*E2*E4*(E5-1) + 0.0014611*E1*E2*(E5-1)^2 + 0.0013839*E1*E4^3 - 0.0012045*E1*E4^2*(E5-1) + 0.00058402*E1*E4*(E5-1)^2 - 0.00025779*E1*(E5-1)^3 - 8.9971E-05*E2^4 - 0.00023374*E2^3*E4 - 0.0004054*E2^3*(E5-1) + 0.0001382*E2^2*E4^2 - 0.00096961*E2^2*E4*(E5-1) - 0.0006973*E2^2*(E5-1)^2 + 0.00036225*E2*E4^3 + 0.000104*E2*E4^2*(E5-1) - 0.00030569*E2*E4*(E5-1)^2 - 0.00010187*E2*(E5-1)^3 - 3.5003E-05*E4^4 - 1.6415E-05*E4^3*(E5-1) + 0.00030839*E4^2*(E5-1)^2 - 0.00026404*E4*(E5-1)^3 + 5.0552E-05*(E5-1)^4 + 0.5685*E1^3 - 0.013487*E1^2*E2 - 0.1014*E1^2*E4 - 0.24111*E1^2*(E5-1) + 0.40244*E1*E2^2 - 0.60635*E1*E2*E4 + 0.096636*E1*E2*(E5-1) + 0.43248*E1*E4^2 + 0.091698*E1*E4*(E5-1) + 0.3174*E1*(E5-1)^2 - 0.11196*E2^3 + 0.079506*E2^2*E4 - 0.13217*E2^2*(E5-1) - 0.041685*E2*E4^2 + 0.057846*E2*E4*(E5-1) - 0.030267*E2*(E5-1)^2 - 0.052105*E4^3 + 0.074374*E4^2*(E5-1) - 0.085138*E4*(E5-1)^2 + 0.046011*(E5-1)^3 + 0.51679*E1^2 - 0.033724*E1*E2 - 0.10772*E1*E4 - 0.12182*E1*(E5-1) + 0.09946*E2^2 + 0.0099693*E2*E4 + 0.030688*E2*(E5-1) + 0.076852*E4^2 + 0.028799*E4*(E5-1) - 0.051907*(E5-1)^2 + 0.3904*E1 - 0.038254*E2 - 0.074777*E4 + 0.0020542*(E5-1) + ubar(1);
    %w(k)=-0.0097589*E1^4 + 0.013093*E1^3*E2 + 0.0076297*E1^3*E4 - 0.013325*E1^3*(E5-1) + 0.0001311*E1^2*E2^2 + 0.0022149*E1^2*E2*E4 + 0.0090055*E1^2*E2*(E5-1) - 0.001577*E1^2*E4^2 - 0.00099654*E1^2*E4*(E5-1) - 0.0042064*E1^2*(E5-1)^2 + 0.001053*E1*E2^3 + 0.0017194*E1*E2^2*E4 - 0.00077488*E1*E2^2*(E5-1) + 0.0030812*E1*E2*E4^2 + 0.0027839*E1*E2*E4*(E5-1) + 0.0037113*E1*E2*(E5-1)^2 + 0.0014337*E1*E4^3 - 0.00090538*E1*E4^2*(E5-1) + 5.1268E-06*E1*E4*(E5-1)^2 - 0.0019199*E1*(E5-1)^3 + 0.0002388*E2^4 + 0.00038528*E2^3*E4 + 0.00056883*E2^3*(E5-1) + 0.00059079*E2^2*E4^2 + 0.00064867*E2^2*E4*(E5-1) + 0.00062928*E2^2*(E5-1)^2 + 0.00055946*E2*E4^3 + 0.0011917*E2*E4^2*(E5-1) + 0.0014356*E2*E4*(E5-1)^2 + 0.0010888*E2*(E5-1)^3 - 6.0349E-05*E4^4 - 0.0003272*E4^3*(E5-1) - 0.00013365*E4^2*(E5-1)^2 - 9.7447E-05*E4*(E5-1)^3 + 0.00021416*(E5-1)^4 + 0.12993*E1^3 + 0.60947*E1^2*E2 - 0.089768*E1^2*E4 - 0.038429*E1^2*(E5-1) - 0.22409*E1*E2^2 + 0.20994*E1*E2*E4 + 0.030654*E1*E2*(E5-1) + 0.018119*E1*E4^2 - 0.2001*E1*E4*(E5-1) + 0.02043*E1*(E5-1)^2 + 0.41838*E2^3 - 0.22373*E2^2*E4 + 0.028704*E2^2*(E5-1) - 0.25451*E2*E4^2 + 0.053264*E2*E4*(E5-1) - 0.23236*E2*(E5-1)^2 - 0.17618*E4^3 + 0.1481*E4^2*(E5-1) - 0.46311*E4*(E5-1)^2 + 0.10457*(E5-1)^3 + 0.056876*E1^2 - 0.15097*E1*E2 - 0.02625*E1*E4 + 0.0078448*E1*(E5-1) + 0.0060111*E2^2 + 0.021899*E2*E4 - 0.06265*E2*(E5-1) - 0.042124*E4^2 + 0.074258*E4*(E5-1) - 0.022069*(E5-1)^2 + 0.086113*E1 + 0.57821*E2 + 0.73259*E4 - 0.064255*(E5-1) + ubar(2);
    
    v(k) = -7.3172E-05*E1^2 - 2.1564E-06*E1*E2 - 1.0276E-05*E1*E4 + 0.00045308*E1*(E5-1) + 0.01617*E2^2 + 0.23249*E2*E4 + 7.412E-07*E2*(E5-1) + 0.0024581*E4^2 - 3.5679E-09*E4*(E5-1) + 0.00019681*(E5-1)^2 + 0.64333*E1 + 5.4303E-05*E2 + 6.2507E-06*E4 + 0.094916*(E5-1) + ubar(1);
    w(k) = 1.0499E-05*E1^2 - 0.11705*E1*E2 - 0.012446*E1*E4 - 7.8688E-07*E1*(E5-1) - 3.052E-06*E2^2 - 1.4085E-06*E2*E4 + 0.01042*E2*(E5-1) + 6.472E-08*E4^2 - 0.12104*E4*(E5-1) - 1.727E-07*(E5-1)^2 - 5.3528E-07*E1 + 0.021046*E2 + 0.19804*E4 + 1.645E-07*(E5-1) + ubar(2);


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