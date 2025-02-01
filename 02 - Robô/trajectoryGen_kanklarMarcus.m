%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mobile robot control on a reference path 
%%% Gregor Klancar, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reprodu��o do artigo, c�digo com ru�do
% 29.12.12

clear all
close all
clc

format long g     


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Inicializa��o de Vari�veis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 0.1;          % sampling time
Tsim = 80;   %90      % tempo de simula��o
h = 0.1;           % integrating time

xa = [0;0;0.0];  % x,y,theta referenciais para in�cio de trajet�ria. 
X1 = xa;         %Euler

xb = [0;0;0.0];
X2 = xb;

% a=1:round(Tsim/h)*0.3;
% b=round(Tsim/h)*0.3:round(Tsim/h)*0.6;
% c=round(Tsim/h)*0.6:round(Tsim/h)*0.8;

% Runge kutta

% Y(n+1) = Y(n) + (1/6)(k1 +2 k2 +2 k3 +k4);
% k1 = h f(X(n),Y(n))
% k2 = h f(X(n)+h/2, Y(n)+k1/2)
% k3 = h f(X(n)+h/2, Y(n) +k2/2)
% k4 = h f(X(n)+h, Y(n)+k3)

% The next loop difines the reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Gera��o da traget�ria e m�todos de seguimento de trajet�ria.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lado = 2;
for k = 1:round(Tsim/h);  % 1 a tempo de simula��o/ tempo de integra��o   
       
       
       if (k<=lado/(Ts*h));   %419
        xb(3)=0;
        xa(3)=0;
     %   x0(3)= 0;
        v(k) = 0.10; %= 0.3;
        w(k) = 0;  % varia��o da velocidade angular positiva
     %  elseif (k>=220)&&(k<=250);
     %   v(k) = 0.00; %= 0.3;
       % w(k) = 0.523; 
      elseif (k>=lado/(Ts*h))&&(k<=2*lado/(Ts*h));
          xb(3) = pi/2;
          xa(3)=pi/2;
          %x0(3)= pi/2;
        v(k) = 0.10; %= 0.3;
        w(k) = 0.0;  % varia��o da velocidade angular positiva
     % elseif (k>=470)&&(k<=500);
      %  v(k) = 0.00; %= 0.3;
      %  w(k) = 0.523;
        elseif (k>=2*lado/(Ts*h))&&(k<=3*lado/(Ts*h));
            xb(3) = pi;
            xa(3)=pi;
          %  x0(3)= pi;
        v(k) = 0.10; %= 0.3;
        w(k) = 0.0;  % varia��o da velocidade angular positiva
     %   elseif (k>=720)&&(k<=750);
      %  v(k) = 0.00; %= 0.3;
      %  w(k) = 0.523;
       elseif (k>=3*lado/(Ts*h))&&(k<=4*lado/(Ts*h));
          xb(3)= 3*pi/2;
          xa(3)=3*pi/2;
        %  x0(3)= 3*pi/2;
        v(k) = 0.10; %= 0.3;
        w(k) = 0.00; % varia��o da velocidade angular negativa 
       else
           v(k)=0;
           w(k)=0;
       end
    
    %Euler
    X1 = [X1, X1(:,k)+ h*([v(k)*cos(xa(3));v(k)*sin(xa(3));w(k)])];  %Euler
    xa = X1(:,k+1);
    %marcus
    X2 = [X2, X2(:,k)+ h*([v(k)*cos(xb(3));v(k)*sin(xb(3));0])];  %Euler
   % xb = X2(:,k+1);
    X2(3,k)=xb(3);
end


Xr = X2;
Ur = [v;w];
%plot(X(1,:),X(2,:),'k--')
hold on
plot(Xr(1,:),Xr(2,:),':')

% grid

vr = v;
wr = w;

% N = 5;
% Q = [1 0 0;
%     0 1 0;
%     0 0 0.1];   %no segundo exemplo foi usado 0.1
% R = 0.1*eye(2);
% vmin = -0.4;
% vmax = 0.4;
% wmin = -0.4;
% wmax = 0.4;
% lb = [];
% ub = [];

% for k = 1:N
%     Qb(3*k-2:3*k,3*k-2:3*k) =Q;
%     Rb(2*k-1:2*k,2*k-1:2*k)=R;
%     lb = [lb;[vmin;wmin]];
%     ub = [ub;[vmax;wmax]]
% end
% 
% % matrizes de pondera��o do Nepsac
% Qn = Qb;
% Rn = [R(1,1)*eye(N),zeros(N);zeros(N),R(2,2)*eye(N)];


%disp('Posi��o inicial')
%disp('[x;y;theta]')
%x0 = input('[x;y;theta] = ')%[1;-1; pi/2];  % x,y,theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Simula��o do rob�.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [-0.5;-0.5;0.0];  % x0 = [0;-1; pi]; cond inicial do robo
X = x0;
v = [];
w = [];

Xm = x0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Coeficientes do controlador 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ksi = 0.3;      %ksi = 0.6; %coeficiente de amortecimento - valor dado no artigo  0.6  Ksi
ohmega_n = 1;   %ohmega_n = 2; %ohmega_n= sqrt(Ur2(t)^2+g*Ur1(t)^2); 
g = 40;          %g = 5;

for k = 1:round(Tsim/h)-1;    % 1 a tempo de simula��o/ tempo de integra��o  
    
    
    a1 = Xr(1,k)-x0(1);
    a2 = Xr(2,k)-x0(2);
    a3 = Xr(3,k)-x0(3);
    
    e1 = cos(x0(3))*a1+sin(x0(3))*a2;
    e2 = -sin(x0(3))*a1+cos(x0(3))*a2;
    e3 = a3;
    
    %controlador K
    
   k1 = 2*ksi*ohmega_n;
   k2 = g*abs(vr(k)); 
   k3 = k1;
%     K = [-k1 0 0; 0 -sign(v(k))*k2 -k3];
%     v = K*e;
%    v1 = v(1,:);
%    v2 = v(2,:);
   
    v1 = -k1*(e1);                     % controlador  
    v2 = - sign(vr(k))*k2*(e2)-k3*e3;  % controlador
    
%     v(k) = -2.4*(e1);
%     w(k) = - (e2)*k2-2.4*e3;

    v(k) = vr(k)*cos(e3) - v1; 
    w(k) = wr(k) -  v2;

    % %restri��es
    % 
    % if v(k) >=0.4  %0.4
    %     v(k)=0.4;
    % end
    % 
    % if v(k)< -0.4;
    %     v(k) = -0.4;
    % end
    % 
    % if w(k)>= 0.8;  %=0.4
    %     w(k) = 0.8;  %=0.4
    % end
    % 
    % if w(k)<-0.8;  %=0.4
    %     w(k)=-0.8;  %=0.4
    % end
    % 
%     %inserir o ru�do 
%     
%     v(k) = r(k)+v(k);
%     w(k) = r(k)+w(k);
    

%      k1 = h*[v(k)*cos(x0(3)); v(k)*sin(x0(3));w(k)];
%      k2 = h*[v(k)*cos(x0(3)+k1(3)/2); v(k)*sin(x0(3)+k1(3)/2);w(k)];
%      k3 = h*[v(k)*cos(x0(3)+k2(3)/2); v(k)*sin(x0(3)+k2(3)/2);w(k)];
%      k4 = h*[v(k)*cos(x0(3)+k3(3)); v(k)*sin(x0(3)+k3(3));w(k)];
%     
%     X = [X, X(:,k) + (1/6)*(k1+2*k2+2*k3+k4)];
    
 %   c�lculo da sa�da do modelo
    
    Xm(1,1) = x0(1) + Ts*v(k)*cos(x0(3));
    Xm(2,1) = x0(2) + Ts*v(k)*sin(x0(3));
    Xm(3,1) = x0(3) + Ts*w(k);
    
    
%     %inserir ruido de medi��o
%    Xm = Xm + r(k);
%     
    X = [X,Xm];
    x0 =Xm;
    
    %x0 = X(:,k+1);
    
    Erro(:,k) = Xr(:,k+1)- X(:,k+1);
%     k
end


plot(X(1,:),X(2,:),'r')
%legend('rungeKutta','euler','Trajet�ria Simulada')
xlabel('x[m]')
ylabel('y[m]')

% print  -depsc trajetoriaruido
% print  -djpeg90 trajetoriaruido
% print  -dpdf trajetoriaruido

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

% print -depsc velocidaderuido
% print -djpeg90 velocidaderuido

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem do erro x y e theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(tempo,Erro(1,:),tempo,Erro(2,:),tempo,Erro(3,:));
legend('x','y','theta')

 E1 = Erro(1,:);
 E2 = Erro(2,:);
 E3 = Erro(3,:);
%ERRO_TOTAL = E1*E1'*Q(1,1)+E2*E2'*Q(2,2)+E3*E3'*Q(3,3);
% grid
% print -depsc erroruido
% print -djpeg90 erroruido
 
%commandwindow  % mostra a pagina principal matlab