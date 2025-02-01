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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Simulacao do robo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [-1.3;0.3;0];
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

    v(k) = -3.263E-15*E1^2 - 0.00020588*E1*E2 - 7.0684E-07*E1*E3 + 0.0084552*E1*E4 - 0.13432*E1*(E5-1) + 2.6886E-05*E2^2 + 0.42981*E2*E3 - 0.029225*E2*E4 - 0.011229*E2*(E5-1) - 3.0579E-14*E3^2 + 0.010532*E3*E4 - 0.00064648*E3*(E5-1) - 0.00022975*E4^2 - 0.017758*E4*(E5-1) + 0.00052467*(E5-1)^2 + 0.96788*E1 - 0.037226*E2 + 0.026957*E3 - 0.20809*E4 + 0.15834*(E5-1) + ubar(1);
    w(k) = 1.0814E-08*E2^2 - 0.7148*E2*E3 + 0.23874*E2*E4 - 7.4749E-14*E3^2 + 0.0097085*E3*E4 - 0.0089937*E4^2 + 0.10883*E2 + 0.037745*E3 + 0.81506*E4 + ubar(2);


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




