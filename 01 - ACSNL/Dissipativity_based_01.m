clc; clear; close all

% % % Recurrent Dissipativity-Based Inequalities for Controller Design
% % Madeira, Diego de S.; Machado, Gabriel.
% Numerical Examples - A: Global Stabilization of a Polynomial System by State Feedback

% Declaração de variáveis
pvar x1 x2 u
vars = [x1;x2;u];
x = [x1;x2];

% Sistema
f = [2*x1^3+x1^2*x2-6*x1*x2^2+5*x2^3;
     0];

g = [0;
     1];

h = sum(monomials(x,3));

y = h;

n = size(f,2);
m = size(g,2);
p = size(h,2);

% Inputs
beta_V = 1e-8;
beta_T = 1e-8;
beta_L = 1e-8;
n_V = 1;
n_T = 2;
n_L = 1;
alpha = 0;   % alpha0

% Número de iterações
k = 1;
kmax = 100;

% STEP 1: Determinar V0, T0, L0, Q0, S0, R0
prog = sosprogram(vars);

[prog, V] = sospolyvar(prog,monomials(x,2:4),'wscoeff');
[prog, T] = sospolyvar(prog,monomials(x,2:4),'wscoeff');
[prog, L] = sospolyvar(prog,monomials(x,0),'wscoeff');
[prog, Q] = sospolymatrixvar(prog,monomials(x,0),[p,p],'symmetric');
[prog, S] = sospolymatrixvar(prog,monomials(x,0),[p,m]);
[prog, R] = sospolymatrixvar(prog,monomials(x,0),[m,m],'symmetric');

gradV = [diff(V,x1);diff(V,x2)];
dissipativity = -gradV'*(f+g*u) - T + h'*Q*h + 2*h'*S*u + u'*R*u - alpha*(1-L);
prog = sosineq(prog,dissipativity);

rad_ilimitada_V = V - beta_V*(x1^2+x2^2)^n_V;
prog = sosineq(prog,rad_ilimitada_V);
rad_ilimitada_T = T - beta_T*(x1^2+x2^2)^n_T;
prog = sosineq(prog,rad_ilimitada_T);
rad_ilimitada_L = L - beta_L*(x1^2+x2^2)^n_L;
prog = sosineq(prog,rad_ilimitada_L);

prog = sosineq(prog,R);


% Solução do STEP 1
sol = sossolve(prog);

Q0 = double(sosgetsol(sol,Q));
R0 = double(sosgetsol(sol,R));
S0 = double(sosgetsol(sol,S));
L = sosgetsol(sol,L);

delta = S0*inv(R0)*S0'-Q0
if delta >= 0
    K = -inv(R0)*S0'
    V = sosgetsol(sol,V)
else
    % STEP 2: Método iterativo para V, T, alpha, Q, S, R
    while k <= kmax
        prog = sosprogram(vars);

        [prog, V] = sospolyvar(prog,monomials(x,2:4),'wscoeff');
        [prog, T] = sospolyvar(prog,monomials(x,2:4),'wscoeff');
        [prog, alpha] = sospolyvar(prog,monomials(x,0),'wscoeff');
        [prog, Q] = sospolymatrixvar(prog,monomials(x,0),[p,p],'symmetric');
        [prog, S] = sospolymatrixvar(prog,monomials(x,0),[p,m]);
        [prog, R] = sospolymatrixvar(prog,monomials(x,0),[m,m],'symmetric');

        gradV = [diff(V,x1);diff(V,x2)];
        dissipativity = -gradV'*(f+g*u) - T + h'*Q*h + 2*h'*S*u + u'*R*u - alpha*(1-L);
        prog = sosineq(prog,dissipativity);

        rad_ilimitada_V = V - beta_V*(x1^2+x2^2)^n_V;
        prog = sosineq(prog,rad_ilimitada_V);
        rad_ilimitada_T = T - beta_T*(x1^2+x2^2)^n_T;
        prog = sosineq(prog,rad_ilimitada_T);
        rad_ilimitada_L = L - beta_L*(x1^2+x2^2)^n_L;
        prog = sosineq(prog,rad_ilimitada_L);
        
        prog = sosineq(prog,R0-R);
        delta_crescente = S*inv(R0)*S0' + S0*inv(R0)*S' - 2*S0*inv(R0)*S0' + Q0 - Q;
        prog = sosineq(prog,delta_crescente);

        prog = sosineq(prog,R);
                        
        % Solução do STEP 2
        sol = sossolve(prog);

        Q = double(sosgetsol(sol,Q))
        R = double(sosgetsol(sol,R))
        S = double(sosgetsol(sol,S))
        
        delta = S*inv(R)*S'-Q
        if delta >= 0 | k == kmax
            K = -inv(R)*S'
            V = sosgetsol(sol,V)
            break
        end
        k=k+1;
        Q0 = Q; R0 = R; S0 = S;
    end
end

a=5;
passo = 0.01;
[x1,x2] = meshgrid(-a:passo:a,-a:passo:a);

dx1 = 2*x1.^3+x1^2.*x2-6.*x1.*x2.^2+5.*x2.^3;
dx2 = 0.*x1;
figure
streamslice(x1,x2,dx1,dx2)
title('Malha aberta'); xlabel('x1'); ylabel('x2');

a=5;
passo = 0.01;
[x1,x2] = meshgrid(-a:passo:a,-a:passo:a);

dx1 = 2*x1.^3+x1^2.*x2-6.*x1.*x2.^2+5.*x2.^3;
dx2 = K.*(x1.^3 + x1.^2.*x2 + x1.*x2.^2 + x2.^3);
figure
streamslice(x1,x2,dx1,dx2)
title('Malha fechada'); xlabel('x1'); ylabel('x2');