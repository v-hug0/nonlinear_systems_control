clc
clear all

x0 = [0.5; 0.5; 0.5];

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % RESPOSTA EM MALHA ABERTA
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t, x] = ode45(@openloop, [0 0.9], x0);
subplot(2,1,1)
response_OL = plot(t,x(:,1),t,x(:,2),t,x(:,3));
title('Malha Aberta')
legend('x1','x2','x3','Location','best')

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % RESPOSTA EM MALHA FECHADA
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t, x] = ode45(@closedloop, [0 20], x0);
subplot(2,1,2)
response_CL = plot(t,x(:,1),t,x(:,2),t,x(:,3));
title('Malha Fechada')
legend('x1','x2','x3','Location','best')


figure
[x1,x2]=meshgrid(-5:0.01:5, -5:0.01:5);
[x1,x3]=meshgrid(-5:0.01:5, -5:0.01:5);
dx1 = -0.878.*x1+x3-x1.^2.*x3-0.0896.*x1.*x3-0.019.*x2.^2+0.473.*x1.^2+3.813.*x1.^3;
dx2 = x2;
dx3 = -4.209.*x1-0.396.*x3-0.408.*x1.^2-2.137.*x1.^3;
subplot(3,1,1)
streamslice(x1,x2,dx1,dx2,4)
subplot(3,1,2)
streamslice(x1,x2,dx1,dx3,4)
subplot(3,1,3)
streamslice(x1,x2,dx2,dx3,4)

figure
[x1,x2]=meshgrid(-0.6:0.01:0.6, -0.6:0.01:0.6);
[x1,x3]=meshgrid(-0.6:0.01:0.6, -0.6:0.01:0.6);
k = - 0.60357.*x1 + 0.51918.*x2 + 2.6414.*x3 + 1.1853.*x1.^2 - 0.0039012.*x2.^2 + 0.15242.*x3.^2 - 0.2542.*x1.*x2 - 0.58084.*x1.*x3 + 0.15914.*x2.*x3 + 8.0966.*x1.^3 + 0.07795.*x2.^3 + 0.25808.*x3.^3 - 1.9417.*x1.^2.*x2 - 4.2947.*x1.^2.*x3 + 0.5879.*x1.*x2.^2 + 0.34048.*x1.*x3.^2 + 0.077199.*x2.^2.*x3 - 0.13242.*x2.*x3.^2 + 0.69357.*x1.*x2.*x3;
dx1 = -0.878.*x1+x3-x1.^2.*x3-0.0896.*x1.*x3-0.019.*x2.^2+0.473.*x1.^2+3.813.*x1.^3 - 0.216*k;
dx2 = x2;
dx3 = -4.209.*x1-0.396.*x3-0.408.*x1.^2-2.137.*x1.^3 - 20.991*k;
%subplot(3,1,1)
streamslice(x1,x2,dx1,dx2,4)
%subplot(3,1,2)
figure
streamslice(x1,x2,dx1,dx3,4)
%subplot(3,1,3)
figure
streamslice(x1,x2,dx2,dx3,4)


function dx = openloop(t,x)
    dx = [  -0.878*x(1)+x(3)-x(1)^2*x(3)-0.0896*x(1)*x(3)-0.019*x(2)^2+0.473*x(1)^2+3.813*x(1)^3;
            x(3);
            -4.209*x(1)-0.396*x(3)-0.408*x(1)^2-2.137*x(1)^3];
end

function dx = closedloop(t,x)
    k = - 0.60357*x(1) + 0.51918*x(2) + 2.6414*x(3) + 1.1853*x(1)^2 - 0.0039012*x(2)^2 + 0.15242*x(3)^2 - 0.2542*x(1)*x(2) - 0.58084*x(1)*x(3) + 0.15914*x(2)*x(3) + 8.0966*x(1)^3 + 0.07795*x(2)^3 + 0.25808*x(3)^3 - 1.9417*x(1)^2*x(2) - 4.2947*x(1)^2*x(3) + 0.5879*x(1)*x(2)^2 + 0.34048*x(1)*x(3)^2 + 0.077199*x(2)^2*x(3) - 0.13242*x(2)*x(3)^2 + 0.69357*x(1)*x(2)*x(3);
    dx = [  -0.878*x(1)+x(3)-x(1)^2*x(3)-0.0896*x(1)*x(3)-0.019*x(2)^2+0.473*x(1)^2+3.813*x(1)^3-0.216*k;
            x(3);
            -4.209*x(1)-0.396*x(3)-0.408*x(1)^2-2.137*x(1)^3-20.991*k];
end

%%

pvar x1 x2 x3 k;
vars = [x1; x2; x3; k];
x = [x1; x2; x3];
prog = sosprogram(vars);

f = [-0.878*x1+x3-x1^2*x3-0.0896*x1*x3-0.019*x2^2+0.473*x1^2+3.813*x1^3;
     x3;
     -4.209*x1-0.396*x3-0.408*x1^2-2.137*x1^3];
g = [-0.216; 0; -20.991];
x = [x1; x2; x3];   % VETOR DE ESTADOS
k = - 0.60357*x1 + 0.51918*x2 + 2.6414*x3 + 1.1853*x1^2 - 0.0039012*x2^2 + 0.15242*x3^2 - 0.2542*x1*x2 - 0.58084*x1*x3 + 0.15914*x2*x3 + 8.0966*x1^3 + 0.07795*x2^3 + 0.25808*x3^3 - 1.9417*x1^2*x2 - 4.2947*x1^2*x3 + 0.5879*x1*x2^2 + 0.34048*x1*x3^2 + 0.077199*x2^2*x3 - 0.13242*x2*x3^2 + 0.69357*x1*x2*x3;

% Monomios   
monDiss = [[monomials(x,[2 3 4 5 6])]]; %graus 2,3,4,5 e 6 em 'x' 
monQ = monomials(x,[1 2 3]);         % monomios da saida h(x)
h = ones(1,length(monQ))*monQ; % h(x), saída de dimensão 1

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERIFICAÇÃO DA QSR-DISSIPATIVIDADE ESTRITA EM MALHA ABERTA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------- RESTRIÇÃO 1: V>0 -----------------------------%
[prog,V] = sospolyvar(prog,[monomials(x,[2 3 4])],'wscoeff'); % V de grau 4
prog = sosineq(prog,V-10^(-1)*(x1^4+x2^4+x3^4));

%-------------------------- RESTRIÇÃO 2: T>0 -----------------------------%
[prog,T] = sospolyvar(prog,[monomials(x,2)],'wscoeff'); % T de grau 2
prog = sosineq(prog,T-10^(-5)*(x1^2+x2^2+x3^2));

%--------------------- Definindo a 'supply rate' -------------------------%
[prog,Q] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(h)]); %Q É UMA MATRIZ CONSTANTE 
[prog,S] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(k)]); %S É UMA MATRIZ CONSTANTE 

%RESTRIÇÃO 3: 1º TERMO DA DISSIPATIVIDADE
grad = [diff(V,x1); diff(V,x2); diff(V,x3)];  % gradiente de V
diss = -(grad'*f + T - h'*Q*h); % condição de grau 6
prog = sosineq(prog,diss-10^-5*(x1^6+x2^6+x3^6)); % DISS com T(x)>0 é SOS
prog = soseq(prog,(0.5*grad'*g - h'*S)); % RESTRIÇÃO DE IGUALDADE
%prog = sossetobj(prog,Q);
prog = sossolve(prog);

% =============================================
% Finally, get solution

V = sosgetsol(prog,V);
SOLV = sosgetsol(prog,V)
[Q1,Z1] = findsos(SOLV)

T = sosgetsol(prog,T);
SOLT = sosgetsol(prog,T)
[QT,ZT] = findsos(SOLT)

Q = sosgetsol(prog,Q)
S = sosgetsol(prog,S)
  

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERIFICAÇÃO DA QSR-DISSIPATIVIDADE ESTRITA EM MALHA FECHADA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prog = sosprogram(vars);

%-------------------------- RESTRIÇÃO 1: V>0 -----------------------------%
[prog,V] = sospolyvar(prog,[monomials(x,[2 3 4])],'wscoeff'); % V de grau 4
prog = sosineq(prog,V-10^(-1)*(x1^4+x2^4+x3^4));

%-------------------------- RESTRIÇÃO 2: T>0 -----------------------------%
[prog,T] = sospolyvar(prog,[monomials(x,2)],'wscoeff'); % T de grau 2
prog = sosineq(prog,T-10^(-5)*(x1^2+x2^2+x3^2));

%--------------------- Definindo a 'supply rate' -------------------------%
[prog,Q] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(h)]); %Q É UMA MATRIZ CONSTANTE 
[prog,S] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(k)]); %S É UMA MATRIZ CONSTANTE 
[prog,R] = sospolymatrixvar(prog,monomials(vars,0),[length(k),length(k)]); %Q É UMA MATRIZ CONSTANTE 

%RESTRIÇÃO 3: 1º TERMO DA DISSIPATIVIDADE
grad = [diff(V,x1); diff(V,x2); diff(V,x3)];  % gradiente de V
diss = -(grad'*(f+g*k) + T - h'*Q*h - 2*h'*S*k - k'*R*k); % condição de grau 6
prog = sosineq(prog,diss-10^-5*(x1^6+x2^6+x3^6)); % DISS com T(x)>0 é SOS

%prog = sossetobj(prog,Q);
prog = sossolve(prog);

% =============================================
% Finally, get solution

V = sosgetsol(prog,V);
SOLV = sosgetsol(prog,V)
[Q1,Z1] = findsos(SOLV)

T = sosgetsol(prog,T);
SOLT = sosgetsol(prog,T)
[QT,ZT] = findsos(SOLT)

Q = sosgetsol(prog,Q)
S = sosgetsol(prog,S)
R = sosgetsol(prog,R)

% % DIAGRAMA DE FASES
% figure
% syms x1 x2 x3
% [x1,x2,x3]=meshgrid(-5:.1:5, -5:.1:5, -5:.1:5);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k = - 0.60357.*x1 + 0.51918.*x2 + 2.6414.*x3 + 1.1853.*x1^2 - 0.0039012.*x2^2 + 0.15242.*x3^2 - 0.2542.*x1.*x2 - 0.58084.*x1.*x3 + 0.15914.*x2.*x3 + 8.0966.*x1^3 + 0.07795.*x2^3 + 0.25808.*x3.^3 - 1.9417.*x1.^2.*x2 - 4.2947.*x1.^2.*x3 + 0.5879.*x1.*x2.^2 + 0.34048.*x1.*x3.^2 + 0.077199.*x2.^2.*x3 - 0.13242.*x2*x3.^2 + 0.69357.*x1.*x2.*x3;
% dx1=-0.878.*x1+x3-x1.^2.*x3-0.0896.*x1.*x3-0.019.*x2.^2+0.473.*x1.^2+3.813.*x1.^3-0.216*k;
% dx2=x3;
% dx3=-4.209.*x1-0.396.*x3-0.408.*x1.^2-2.137.*x1.^3-20.991*k;
% streamslice(x1,x2,x3,dx1,dx2,dx3);
% hold on
% plot(0,0,'o')
% hold on
% plot(0,0,'x')
% title('Phase portrait')
