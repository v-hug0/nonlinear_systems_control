clc;
clear;
close all;

pvar x1 x2 u;
var = [x1;x2;u];

f = [2*x1^3+x1^2*x2-6*x1*x2^2+5*x2^3; 0];
g = [0;1];
h = x1^3+x1^2*x2+x1*x2^2+x2^3;

n = 2; %2 estados
p = 1; %numero de colunas do h
m = 1; %numero de colunas do g

beta = 1e-8;
nv = 1;
nt = 2;


k = 0;
kmax = 100; 

prog=sosprogram(var);

[prog,V0] = sospolyvar(prog,monomials([x1;x2],2:6),'wscoeff');
[prog,T0] = sospolyvar(prog,monomials([x1;x2],2:6),'wscoeff');
[prog,Q0]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
[prog,S0]=sospolymatrixvar(prog,monomials(var,0),[p,m]); 
[prog,R0]=sospolymatrixvar(prog,monomials(var,0),[m,m],'symmetric');

gradV0 = [diff(V0,x1); diff(V0,x2)];
cond1 = -gradV0'*(f+g*u)-T0+h'*Q0*h+2*h'*S0*u+u'*R0*u;
cond2 = V0-beta*(x1^2+x2^2)^nv;
cond3 = T0-beta*(x1^2+x2^2)^nt;
prog = sosineq(prog, cond1);
prog = sosineq(prog, cond2);
prog = sosineq(prog, cond3);
prog = sosineq(prog, R0);

sol = sossolve(prog);

Q0 = double(sosgetsol(sol,Q0));
R0 = double(sosgetsol(sol,R0));
S0 = double(sosgetsol(sol,S0));

if S0*inv(R0)*S0'-Q0 >=0
    K = -inv(R0)*S0';
    V0 = sosgetsol(sol,V0)
else while k<kmax
    prog=sosprogram(var);
    [prog,V] = sospolyvar(prog,monomials([x1;x2],2:6),'wscoeff');
    [prog,T] = sospolyvar(prog,monomials([x1;x2],2:6),'wscoeff');
    [prog,Q]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
    [prog,S]=sospolymatrixvar(prog,monomials(var,0),[p,m]); 
    [prog,R]=sospolymatrixvar(prog,monomials(var,0),[m,m],'symmetric');

    gradV = [diff(V,x1); diff(V,x2)];
    cond1 = -gradV'*(f+g*u)-T+h'*Q*h+2*h'*S*u+u'*R*u;
    cond2 = V-beta*(x1^2+x2^2)^nv;
    cond3 = T-beta*(x1^2+x2^2)^nt;
    cond4 = R0-R;
    cond5 = S*inv(R0)*S0'+S0*inv(R0)*S'-2*S0*inv(R0)*S0'+Q0-Q;
    prog = sosineq(prog, cond1);
    prog = sosineq(prog, cond2);
    prog = sosineq(prog, cond3);
    prog = sosineq(prog, R);
    prog = sosineq(prog, cond4);
    prog = sosineq(prog, cond5);
    %prog = sosineq(prog, cond6);

    sol = sossolve(prog);

    Q = double(sosgetsol(sol,Q));
    R = double(sosgetsol(sol,R));
    S = double(sosgetsol(sol,S));

    Q0 = Q;
    S0 = S;
    R0 = R;

    k=k+1;
    if S*inv(R)*S'-Q >= 0
        K = -inv(R)*S';
        V = sosgetsol(sol,V)
        break
    end
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
