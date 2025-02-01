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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Coeficientes do controlador 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pvar e1 e2 e3 e4 e5
e_n=[e1;e2;e3;e4;e5];
Ur = [0.2 0.2];
f = [Ur(1)*(e5+1); Ur(1)*e4; Ur(2); Ur(2)*(e5+1); -Ur(2)*e4];
g = [-1 e2; 0 -e1; 0 -1; 0 -e5-1; 0 e4];

n = size(f,2);
m = size(g,2);

% V inicial 
V = e_n'*e_n;
ubar = [0.2;0.2];
i = 0;
i_max = 20;

while i<i_max
    % solving with SOSTOOLS (V fixed)
    prog = sosprogram(e_n);
    [prog,kx1] = sospolyvar(prog,[monomials(e_n,1:2)],'wscoeff');
    [prog,kx2] = sospolyvar(prog,[monomials([e2;e3;e4],1:2)],'wscoeff');
    kx=[kx1;kx2];
    [prog,T] = sospolyvar(prog,[monomials(e_n,1:3)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,1:3)],'wscoeff');
    g1 = e4^2+e5^2-1;
    % SOS constraints
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*(kx))-T-s1*g1;
    prog = sosineq(prog,dV);
    % first solution (k and lambda)
    sol = sossolve(prog);
    kx = sosgetsol(sol,kx);
    % solving with SOSTOOLS (k fixed)
    prog = sosprogram(e_n);
    [prog,V] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    [prog,T] = sospolyvar(prog,[monomials(e_n,1:3)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,1:3)], 'wscoeff');
    g1 = e4^2+e5^2-1;
    % SOS constraints
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*(kx))-T-s1*g1;
    prog = sosineq(prog,dV);
    % first solution (V)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    i = i+1
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Simulacao do robo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [-1.5;0.5;pi/10];
x = xo;
X = x;
v = [];
w = [];
ubar = [0.2;0.2];

%limitV = 0.3;
%limitW = 0.6;

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
    
    %o melhor (OLD)
    %v(k) = 2.2178E-14*E1^2 - 3.1002E-06*E1*E2 + 0.00018347*E1*E3 - 1.3284E-05*E1*E4 + 0.00012021*E1*(E5-1) + 8.8884E-06*E2^2 - 3.8116E-05*E2*E3 - 0.00015199*E2*E4 + 0.00045948*E2*(E5-1) - 4.6418E-07*E3^2 + 0.00048292*E3*E4 - 0.00014007*E3*(E5-1) + 0.062241*E4^2 + 0.0017243*E4*(E5-1) + 0.061997*(E5-1)^2 + 1.9108*E1 - 0.028616*E2 - 0.010646*E3 - 0.01574*E4 + 0.15933*(E5-1) + ubar(1);
    %w(k) = -0.0012161*E2^2 + 0.0024543*E2*E3 + 0.0028334*E2*E4 - 9.0476E-05*E3^2 - 0.06354*E3*E4 + 0.00051793*E4^2 + 0.42684*E2 + 0.76584*E3 - 0.02847*E4 + ubar(2);

    % dissipatividade
    %v(k) = -0.0051738*E1^2 - 0.0021226*E1*E2 - 0.026702*E1*E3 + 0.2642*E1*E4 + 1.2022*E1*(E5-1) - 0.018476*E2^2 + 0.0083082*E2*E3 + 0.04823*E2*E4 + 0.016441*E2*(E5-1) - 0.0040279*E3^2 + 0.0065559*E3*E4 - 0.0061168*E3*(E5-1) - 0.023727*E4^2 - 0.013772*E4*(E5-1) - 0.046528*(E5-1)^2 + 4.0633*E1 - 0.50829*E2 + 0.1346*E3 + 0.08552*E4 - 0.060523*(E5-1) + ubar(1);
    %w(k) = 0.91726*E2^2 - 0.13194*E2*E3 - 0.070694*E2*E4 - 0.044194*E3^2 + 0.021534*E3*E4 + 0.016843*E4^2 - 0.016826*E2 + 3.5761*E3 + 0.91726*E4 + ubar(2);

    v(k) = - 0.055407*E1*E4 + 0.024284*E1*(E5-1) + 9.2462E-07*E2^2 - 0.0011518*E2*E3 + 0.019062*E2*E4 - 0.0043647*E2*(E5-1) - 9.5563E-08*E3^2 - 0.02301*E3*E4 + 0.00079743*E3*(E5-1) - 0.16727*E4^2 + 0.0027039*E4*(E5-1) - 0.36928*(E5-1)^2 + 3.4023*E1 - 0.35822*E2 - 0.83481*E3 + 0.10252*E4 + 0.0020579*(E5-1) + ubar(1);
    w(k) = - 0.015316*E2*E4 + 0.017616*E3*E4 - 0.25847*E4^2 + 0.27795*E2 + 0.63668*E3 - 0.041518*E4 + ubar(2);

    % if v(k) >= limitV
    %      v(k) = limitV;
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




