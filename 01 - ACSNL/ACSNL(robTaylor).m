clc
clear all

e0 = [-0.402; -0.333; 0.391];
vr = 0.2;
wr = 0.2;
k1 = 2.6;
k2 = 1;
k3 = 3;

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % RESPOSTA EM MALHA ABERTA
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[e1,e2,e3]=meshgrid(-1:0.5:1, -1:0.5:1, -1:0.5:1);
d_e1 = vr*(1-0.5.*e3.^2+0.042.*e3.^4-0.0014.*e3.^6);
d_e2 = vr*(e3-0.17.*e3.^3+0.0083.*e3.^5-0.0002.*e3.^7);
d_e3 = wr.*ones(size(e3));
quiver3(e1,e2,e3,d_e1,d_e2,d_e3)
hold on
e0 = zeros(3,15);
for i = 1:15
    e0(:,i)=2*rand(3,1)-[1;1;1];
    [t, e] = ode45(@openloop, [0 8], e0(:,i));
    response_CL = plot3(e(:,1),e(:,2),e(:,3));
end
response_OL = plot3(e(:,1),e(:,2),e(:,3));
title('Malha Aberta')

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % RESPOSTA EM MALHA FECHADA
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
[e1,e2,e3]=meshgrid(-1:0.5:1, -1:0.5:1, -1:0.5:1);
K1 = vr*(1-0.5.*e3.^2+0.042.*e3.^4-0.0014.*e3.^6)-k3.*e3.*(wr+0.5*vr*(k2.*(e2+k3.*e3)+(1/k3).*(e3-0.17.*e3.^3+0.0083.*e3.^5-0.0002.*e3.^7)))+k1.*e1;
K2 = wr+0.5*vr*(k2.*(e2+k3.*e3)+(1/k3).*(e3-0.17.*e3.^3+0.0083.*e3.^5-0.0002.*e3.^7));
d_e1 = vr*(1-0.5.*e3.^2+0.042.*e3.^4-0.0014.*e3.^6)-K1+e(2).*K2;
d_e2 = vr*(e3-0.17.*e3.^3+0.0083.*e3.^5-0.0002.*e3.^7)-e(1).*K2;
d_e3 = (wr-K2);
quiver3(e1,e2,e3,d_e1,d_e2,d_e3)
hold on
%e0 = [-1; -1; 0.4];
e0 = zeros(3,15);
for i = 1:15
    e0(:,i)=2*rand(3,1)-[1;1;1];
    [t, e] = ode45(@closedloop, [0 50], e0(:,i));
    response_CL = plot3(e(:,1),e(:,2),e(:,3));
end
title('Malha Fechada')

function de = openloop(t,e)
    vr = 0.2;
    wr = 0.2;
    de = [  vr*(1-0.5*e(3)^2+0.042*e(3)^4-0.0014*e(3)^6);
            vr*(e(3)-0.17*e(3)^3+0.0083*e(3)^5-0.0002*e(3)^7);
            wr];
end

function de = closedloop(t,e)
    vr = 0.2;
    wr = 0.2;
    k1 = 2.6;
    k2 = 1;
    k3 = 3;
    k = [vr*(1-0.5*e(3)^2+0.042*e(3)^4-0.0014*e(3)^6)-k3*e(3)*(wr+0.5*vr*(k2*(e(2)+k3*e(3))+(1/k3)*(e(3)-0.17*e(3)^3+0.0083*e(3)^5-0.0002*e(3)^7)))+k1*e(1);
         wr+0.5*vr*(k2*(e(2)+k3*e(3))+(1/k3)*(e(3)-0.17*e(3)^3+0.0083*e(3)^5-0.0002*e(3)^7))];
    de = [  vr*(1-0.5*e(3)^2+0.042*e(3)^4-0.0014*e(3)^6)-k(1)+e(2)*k(2);
            vr*(e(3)-0.17*e(3)^3+0.0083*e(3)^5-0.0002*e(3)^7)-e(1)*k(2);
            wr-k(2)];
end

%%

vr = 0.2;
wr = 0.2;
k1 = 2.6;
k2 = 1;
k3 = 3;

pvar e1 e2 e3 u;
vars = [e1; e2; e3; u];
e = [e1; e2; e3];
prog = sosprogram(vars);

f = [vr*(1-0.5*e3^2+0.042*e3^4-0.0014*e3^6);
     vr*(e3-0.17*e3^3+0.0083*e3^5-0.0002*e3^7);
     wr];
g = [-1 e2; 0 -e1; 0 -1];
u = [vr*(1-0.5*e3^2+0.042*e3^4-0.0014*e3^6)-k3*e3*(wr+0.5*vr*(k2*(e2+k3*e3)+(1/k3)*(e3-0.17*e3^3+0.0083*e3^5-0.0002*e3^7)))+k1*e1;
     wr+0.5*vr*(k2*(e2+k3*e3)+(1/k3)*(e3-0.17*e3^3+0.0083*e3^5-0.0002*e3^7))];

% Monomios   
monDiss = [[monomials(e,[2 3 4 5 6])]]; %graus 2,3,4,5 e 6 em 'x' 
monQ = monomials(e,[1 2 3]);         % monomios da saida h(x)
h = ones(1,length(monQ))*monQ; % h(x), saída de dimensão 1

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERIFICAÇÃO DA QSR-DISSIPATIVIDADE ESTRITA EM MALHA ABERTA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------- RESTRIÇÃO 1: V>0 -----------------------------%
[prog,V] = sospolyvar(prog,[monomials(e,[2 3 4])],'wscoeff'); % V de grau 4
prog = sosineq(prog,V-10^(-1)*(e1^4+e2^4+e3^4));

%-------------------------- RESTRIÇÃO 2: T>0 -----------------------------%
[prog,T] = sospolyvar(prog,[monomials(e,2)],'wscoeff'); % T de grau 2
prog = sosineq(prog,T-10^(-5)*(e1^2+e2^2+e3^2));

%--------------------- Definindo a 'supply rate' -------------------------%
[prog,Q] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(h)]); %Q É UMA MATRIZ CONSTANTE 
[prog,S] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(u)]); %S É UMA MATRIZ CONSTANTE 

%RESTRIÇÃO 3: 1º TERMO DA DISSIPATIVIDADE
grad = [diff(V,e1); diff(V,e2); diff(V,e3)];  % gradiente de V
diss = -(grad'*f + T - h'*Q*h); % condição de grau 6
prog = sosineq(prog,diss-10^-5*(e1^6+e2^6+e3^6)); % DISS com T(x)>0 é SOS
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
[prog,V] = sospolyvar(prog,[monomials(e,[2 3 4])],'wscoeff'); % V de grau 4
prog = sosineq(prog,V-10^(-1)*(e1^4+e2^4+e3^4));

%-------------------------- RESTRIÇÃO 2: T>0 -----------------------------%
[prog,T] = sospolyvar(prog,[monomials(e,2)],'wscoeff'); % T de grau 2
prog = sosineq(prog,T-10^(-5)*(e1^2+e2^2+e3^2));

%--------------------- Definindo a 'supply rate' -------------------------%
[prog,Q] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(h)]); %Q É UMA MATRIZ CONSTANTE 
[prog,S] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(u)]); %S É UMA MATRIZ CONSTANTE 
[prog,R] = sospolymatrixvar(prog,monomials(vars,0),[length(u),length(u)]); %Q É UMA MATRIZ CONSTANTE 

%RESTRIÇÃO 3: 1º TERMO DA DISSIPATIVIDADE
grad = [diff(V,e1); diff(V,e2); diff(V,e3)];  % gradiente de V
diss = -(grad'*(f+g*u) + T - h'*Q*h - 2*h'*S*u - u'*R*u); % condição de grau 6
prog = sosineq(prog,diss-10^-5*(e1^6+e2^6+e3^6)); % DISS com T(x)>0 é SOS

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
