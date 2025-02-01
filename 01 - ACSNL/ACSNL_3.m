

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            EXEMPLO MOEZ 2009                            %
%                        DISSIPATIVIDADE DELTA>=0                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

syms x1 x2 u real;
vars = [x1; x2; u];

prog = sosprogram(vars);

f = [-x1+x2+x1^2+x1*x2-x1^3+x1^2*x2-x1*x2^2+2*x2^3;
     -x1+1.5*x2-x1^2-0.5*x1*x2-x1^3-x1^2*x2+0.5*x1*x2^2-2*x2^3];
g = [0; 1];
x = [x1; x2];   % VETOR DE ESTADOS

% Monomios   
monDiss = [[monomials([x1; x2],[2 3 4 5 6])]]; %graus 2,3,4,5 e 6 em 'x' 
monQ = monomials([x1; x2],[1 2 3]);         % monomios da saida h(x)
h = ones(1,length(monQ))*monQ; % h(x), saída de dimensão 1


%-------------------------- RESTRIÇÃO 1: V>0 -----------------------------%
[prog,V] = sospolyvar(prog,[monomials([x1; x2],[2 3 4])],'wscoeff'); % V de grau 4
prog = sosineq(prog,V-10^-1*(x1^4+x2^4));

%-------------------------- RESTRIÇÃO 2: T>0 -----------------------------%
[prog,T] = sospolyvar(prog,[monomials([x1; x2],[2])],'wscoeff'); % T de grau 2
prog = sosineq(prog,T-10^-5*(x1^2+x2^2));

%--------------------- Definindo a 'supply rate' -------------------------%
[prog,Q] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(h)]); %Q É UMA MATRIZ CONSTANTE 
[prog,S] = sospolymatrixvar(prog,monomials(vars,0),[length(h),length(u)]); %S É UMA MATRIZ CONSTANTE 

%RESTRIÇÃO 3: 1º TERMO DA DISSIPATIVIDADE
[prog,diss] = sospolyvar(prog,[monDiss],'wscoeff');
grad = [diff(V,x1); diff(V,x2)];  % gradiente de V
diss = -(grad'*f + T - h'*Q*h); % condição de grau 6
prog = sosineq(prog,diss-10^-5*(x1^6+x2^6)); % DISS com T(x)>0 é SOS

prog = soseq(prog,(0.5*grad'*g - h'*S)); % RESTRIÇÃO DE IGUALDADE
 
% prog = sosmatrixineq(prog,Q2); % Q2>0

HessV = [diff(grad(1),x1) diff(grad(1),x2); diff(grad(2),x1) diff(grad(2),x2)];
prog = sosineq(prog,x'*HessV*x-10^-5*(x1^2+x2^2)); % Hessina de V é SOS -> V é fortemente convexa

prog = sossetobj(prog,Q);

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

grad = [diff(V,x1); diff(V,x2)];
diss = -(grad'*f + T - h'*Q*h);
SOLD = sosgetsol(prog,diss)
[Q3,Z3] = findsos(SOLD)

SOLHESSV = sosgetsol(prog,HessV)
[Q4,Z4] = findsos(x'*SOLHESSV*x-10^-5*(x1^2+x2^2)) 



% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % DIAGRAMA DE FASES
% syms x1 x2
% % [x1,x2]=meshgrid(-1.2:.1:1.2, -1.2:.1:1.2);
% % [x1,x2]=meshgrid(0:.01:1.2, 0:.01:0.4);
% [x1,x2]=meshgrid(-20:1:20, -20:1:20);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K = -7.1265;
% u = K*(x1.^3 + (x1.^2).*x2 + x1.^2 + x1.*(x2.^2) + x1.*x2 + x1 + x2.^3 + x2.^2 + x2);
% dx1=-x1+x2+x1.^2+x1.*x2-x1.^3+(x1.^2).*x2-x1.*(x2.^2)+2*x2.^3;
% dx2=-x1+1.5*x2-x1.^2-0.5*x1.*x2-x1.^3-(x1.^2).*x2+0.5*x1.*(x2.^2)-2*x2.^3 + u;
% streamslice(x1,x2,dx1,dx2);
% hold on
% plot(0,0,'o')
% hold on
% plot(0,0,'x')
% title('Phase portrait')
% axis tight equal








 








