% Tarefa 02: Projeto de controlador
%
% Autor: João Mateus. 
%
% a) QSR - Dissipatividade
% b) V-K, variação 01;
% c) V-K, variação 02.
%

%% Início e definição das constantes.
clc
clear 
close all

% Variáveis para os estados
pvar x1 x2 x1n x2n;

% Constantes do motor
L = 6.008;           % La + Lf
R = 0.12;            % Ra + Rf
B = 0.0002;
J = 0.05;
Laf = 1.766;
% V = 0 ~ 48;
i = 27.852447355566078;         % 78; i(max) = 250
w = 2800;                       % P/ V = 25 km/h

% Constantes do veículo 
m = 800;
A = 1.8;
p = 1.25;
Cd = 0.3;
phi = 0;
urr = 0.015;
r = 0.25;
G = 11;
grav = 9.88;

% Constantes
a1 = R/L;
a2 = Laf/L;
a3 = Laf/(J+m*(r^2/G^2));
a4 = B/(J+m*(r^2/G^2));
a5 = (r*urr*m*grav)/(G*(J+m*(r^2/G^2)));
a6 = -(p*A*Cd*r^2)/(2*G^2*(J+m*(r^2/G^2)));
a7 = 1/L;
a8 = -1/((J+m*(r^2/G^2))^2);
a9 = (p*A*Cd*r^2)/G^2;  
a10 = ((r*urr*m*grav)/G);
a11 = 0.5*p*A*Cd*(r^2/G^2);
a12 = J+m*r^2/G^2;


%% Sistema não linear
% Estados
X = [x1; x2];

%Funções do sistema
finicial = [ -a1*x1-a2*x1*x2;
             a3*x1^2-a4*x2-a5+a6*x2^2];
g = [ a7; 0];
h = x2;

% Transladando para a origem
ubarra = 1.377281239774858e+05;

f = [ -a1*(x1+i)-a2*(x1+i)*(x2+w); a3*(x1+i)^2-a4*(x2+w)-a5+a6*(x2+w)^2]+g*ubarra;

% Condições iniciais
V = X'*X;
rho = 1;

%% Inicio do Loop:
% Variáveis auxiliares
n = 0;
feasu = 0;
feasv = 0;

% Loop de interações
while n < 100 

    % Definindo U e alpha;
    prog = sosprogram(X);
    
    [prog, u] = sospolyvar(prog, monomials(X, 1:3), 'wscoeff');

    [prog, alpha] = sospolyvar(prog, monomials(X, 2:4), 'wscoeff');
    prog = sosineq(prog,alpha-1e-6*(x1^2+x2^2));

    gradV = [diff(V, x1); diff(V, x2)];
    cond = -gradV'*(f+g*u)-alpha*(rho-V);
    prog = sosineq(prog, cond-1e-6*(x1^2+x2^2));

    %Pegando a solução;
    prog = sossolve(prog);
    u = sosgetsol(prog,u)
    alpha = sosgetsol(prog,alpha)
    feasu = prog.solinfo.info.feasratio;

    % Encontrando V e rho;
    prog = sosprogram(X);

    [prog, V] = sospolyvar(prog, monomials(X, 2:4), 'wscoeff');
    prog = sosineq(prog,V-1e-6*(x1^2+x2^2));

    [prog, rho] = sospolyvar(prog, monomials(X, 0), 'wscoeff');
    prog = sosineq(prog,rho-1e-6);

    gradV = [diff(V, x1); diff(V, x2)];
    cond = -gradV'*(f+g*u)-alpha*(rho-V);
    prog = sosineq(prog, cond-1e-4*(x1^2+x2^2));

    % Pegando a solução;
    prog = sossolve(prog);
    V = sosgetsol(prog,V)
    rho = sosgetsol(prog,rho)
    feasv = prog.solinfo.info.feasratio;

    % Atualizando o while
    n = n + 1;

    % Feas perto de 1
    if abs(1-feasu) < 0.001 && abs(1-feasv) < 0.001
        break;
    end
end