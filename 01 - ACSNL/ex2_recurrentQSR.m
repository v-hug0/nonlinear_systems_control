clc;
clear;
close all;

pvar x1 x2 u;
var = [x1;x2;u];

f = [-x1+x1*x2;x1+2*x2+x1^2+x1^2*x2];
g = [x2;1];
h = x2+x2^2+x2^3;

n = 2; %2 estados
p = 1; %numero de colunas do h
m = 1; %numero de colunas do g

betaV = 1e-4;
betaT = 1e-4;
betaL = 1e-6;

ni = 1;

alfa0 = 1e-6*(x1^2+x2^2);

k = 0;
kmax = 300; 

prog=sosprogram(var);

[prog,V0] = sospolyvar(prog,monomials([x1;x2],2),'wscoeff');
[prog,T0] = sospolyvar(prog,monomials([x1;x2],2:4),'wscoeff');
[prog,L0] = sospolyvar(prog,monomials([x1;x2],2:4),'wscoeff');
[prog,Q0]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
[prog,S0]=sospolymatrixvar(prog,monomials(var,0),[p,m]); 
[prog,R0]=sospolymatrixvar(prog,monomials(var,0),[m,m],'symmetric');

gradV0 = [diff(V0,x1); diff(V0,x2)];
cond1 = -gradV0'*(f+g*u)-T0+h'*Q0*h+2*h'*S0*u+u'*R0*u-alfa0*(1-L0);
cond2 = V0-betaV*(x1^2+x2^2)^ni;
cond3 = T0-betaT*(x1^2+x2^2)^ni;
cond4 = L0-betaL*(x1^2+x2^2)^ni;
prog = sosineq(prog, cond1);
prog = sosineq(prog, cond2);
prog = sosineq(prog, cond3);
prog = sosineq(prog, cond4);
prog = sosineq(prog, R0);

sol = sossolve(prog);

Q0 = double(sosgetsol(sol,Q0));
R0 = double(sosgetsol(sol,R0));
S0 = double(sosgetsol(sol,S0));
L0 = sosgetsol(sol,L0);

if S0*inv(R0)*S0'-Q0 >=0
    K = -inv(R0)*S0';
    V0 = sosgetsol(sol,V0);
else while k<kmax
    prog=sosprogram(var);
    [prog,V] = sospolyvar(prog,monomials([x1;x2],2),'wscoeff');
    [prog,alfa] = sospolyvar(prog,monomials([x1;x2],2:8),'wscoeff');
    [prog,T] = sospolyvar(prog,monomials([x1;x2],2:4),'wscoeff');
    [prog,Q]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
    [prog,S]=sospolymatrixvar(prog,monomials(var,0),[p,m]); 
    [prog,R]=sospolymatrixvar(prog,monomials(var,0),[m,m],'symmetric');

    gradV = [diff(V,x1); diff(V,x2)];
    cond1 = -gradV'*(f+g*u)-T+h'*Q*h+2*h'*S*u+u'*R*u-alfa*(1-L0);
    cond2 = V-betaV*(x1^2+x2^2)^ni;
    cond3 = T-betaT*(x1^2+x2^2)^ni;
    cond4 = L0-betaL*(x1^2+x2^2)^ni;
    cond5 = R0-R;
    cond6 = S*inv(R0)*S0'+S0*inv(R0)*S'-2*S0*inv(R0)*S0'+Q0-Q;
    prog = sosineq(prog, cond1);
    prog = sosineq(prog, cond2);
    prog = sosineq(prog, cond3);
    prog = sosineq(prog, R);
    prog = sosineq(prog, cond4);
    prog = sosineq(prog, cond5);
    prog = sosineq(prog, cond6);


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
        V = sosgetsol(sol,V);
        break
    end
end
end

%% 
ro = 1e-8;
while true
    prog=sosprogram(var);
    [prog,gamma] = sospolyvar(prog,monomials([x1;x2],2:8),'wscoeff');
    cond = 1-L0-gamma*(ro-V);
    prog = sosineq(prog,cond);
    sol = sossolve(prog);
    gamma = sosgetsol(sol,gamma);
    if abs(1-sol.solinfo.info.feasratio)>=0.2
        ro=ro-1e-5;    
        break
    end
    ro=ro+1e-5;
end

%% 

a=1;
passo = 0.01;
[x1,x2] = meshgrid(-a:passo:a,-a:passo:a);
% malha aberta
dx1 = -x1+x1.*x2;
dx2 = x1+2.*x2+x1.^2+x1.^2.*x2;
figure
streamslice(x1,x2,dx1,dx2)
title('Malha aberta'); xlabel('x1'); ylabel('x2');
figure
%contorno
% V_ = 0.0015177.*x1.^2 - 0.0017249.*x1.*x2 + 0.0029903.*x2.^2;
% L0_= 3.4245.*x1.^4 - 0.0086682.*x1.^3.*x2 + 0.40585.*x1.^2.*x2.^2 + 0.0004442.*x1.*x2.^3 + 0.84634.*x2.^4 - 0.0073863.*x1.^3 + 0.010035.*x1.^2.*x2 + 0.00028607.*x1.*x2.^2 - 0.00042428.*x2.^3 + 0.82866.*x1.^2 - 2.3834e-05.*x1.*x2 + 0.82551.*x2.^2;
% figure
% hold on;
% contour(x1,x2,V_,[ro ro]);
% contour(x1,x2,L0_,[1 1],"EdgeColor",'red');
% title('Malha aberta'); xlabel('x1'); ylabel('x2');
% malha fechada
dx1 = -x1+x1.*x2+x2.*K.*(x2+x2.^2+x2.^3);
dx2 = x1+2.*x2+x1.^2+x1.^2.*x2+K.*(x2+x2.^2+x2.^3);
streamslice(x1,x2,dx1,dx2,2)
title('Malha fechada'); xlabel('x1'); ylabel('x2');
hold on;

elipV = matlabFunction(p2s(ro-V));
zhandle = fimplicit(elipV)
zhandle.LineWidth = 1.5;

elipL = matlabFunction(p2s(1-L0));
zhandle = fimplicit(elipL)
zhandle.LineWidth = 1.5;
zhandle.Color = "b";

