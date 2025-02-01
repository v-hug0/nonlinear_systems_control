%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mobile robot control on a reference path 
%%% Gregor Klancar, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modificado por Marcus
% 08.04.2016

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Coeficientes do controlador 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pvar e1 e2 e3 e4 e5 u1 u2
dpvar S1 S2 S3 S4 S5 S6 S7 S8 S9
var = [e1;e2;e3;e4;e5;u1;u2];
e_n = [e1;e2;e3;e4;e5];
u = [u1;u2];
Ur = [0.2 0.2];
f = [Ur(1)*(e5+1); Ur(1)*e4; Ur(2); Ur(2)*(e5+1); -Ur(2)*e4];
g = [-1 e2; 0 -e1; 0 -1; 0 -e5-1; 0 e4];
h = monomials(e_n,1:2);
%h2 = monomials([e2;e3;e4],1:2);
%h = [h1;h2];

n = size(f,2);
m = size(g,2);
p = size(h,1); %numero de linhas do h
j = size(h,2);

beta = 1e-6;
nv = 2;
nt = 2;

i = 0;
imax = 100; 

prog=sosprogram(var);
[prog,V0] = sospolyvar(prog,monomials(e_n,1:4),'wscoeff');
[prog,T0] = sospolyvar(prog,monomials(e_n,1:4),'wscoeff');
[prog,s1] = sospolyvar(prog, [monomials(e_n,1:4)],'wscoeff');
[prog,Q0]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
[prog,S10]=sospolymatrixvar(prog,monomials(var,0),[p,1]); 
[prog,S20]=sospolymatrixvar(prog,monomials(var,0),[9,1]);


S0 = [S10 [0; S20(1); S20(2); S20(3); 0; 0; 0; 0; 0; 0; 0; 0; 0; S20(3); 0; 0; 0; 0; 0; 0]];
[prog,R10]=sospolymatrixvar(prog,monomials(var,0),[m,1]);
R0 = [R10(1) 0; 0 R10(2)];
g1 = e4^2+e5^2-1;

gradV0 = [diff(V0,e1) diff(V0,e2) diff(V0,e3) diff(V0,e4) diff(V0,e5)];
cond1 = -gradV0*(f+g*u)-T0-s1*g1+h'*Q0*h+2*h'*S0*u+u'*R0*u;
%cond2 = V0-beta*(e1^2+e2^2+e3^2+e4^2+e5^2)^nv;
%cond3 = T0-beta*(e1^2+e2^2+e3^2+e4^2+e5^2)^nt;
%prog = sosineq(prog, cond1);
%prog = sosineq(prog, cond2);
%prog = sosineq(prog, cond3);

sol = sossolve(prog);

Q0 = double(sosgetsol(sol,Q0));
R0 = double(sosgetsol(sol,R0));
S0 = double(sosgetsol(sol,S0));

if S0*inv(R0)*S0'-Q0 >=0
    K = -inv(R0)*S0';
    V0 = sosgetsol(sol,V0);
else 
    while i<imax
    prog=sosprogram(var);
    [prog,V] = sospolyvar(prog,monomials(e_n,1:4),'wscoeff');
    [prog,T] = sospolyvar(prog,monomials(e_n,1:4),'wscoeff');
    [prog,s1] = sospolyvar(prog,monomials(e_n,1:4),'wscoeff');
    [prog,Q]=sospolymatrixvar(prog,monomials(var,0),[p,p],'symmetric');
    [prog,S1]=sospolymatrixvar(prog,monomials(var,0),[p,1]); 
    [prog,S2]=sospolymatrixvar(prog,monomials(var,0),[9,1]);
    S = [S1 [0; S2(1); S2(2); S2(3); 0; 0; 0; S2(3); 0; S2(5); S2(6); 0; S2(7); S2(8); S2(9); 0; 0; 0; 0; 0]];
    %S = [[S1(1); S1(2); S1(3); S1(4); S1(5); 0; 0; 0; 0; 0; 0; 0; 0; 0; S1(6); 0; 0; 0; 0; S1(7)] [0; S2(1); S2(2); S2(3); 0; 0; 0; 0; 0; 0; 0; 0; 0; S2(3); 0; 0; 0; 0; 0; 0]];
    [prog,R1]=sospolymatrixvar(prog,monomials(var,0),[m,1]); 
    R = [R1(1) 0; 0 R1(2)];
    g1 = e4^2+e5^2-1;
    
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    cond1 = -gradV*(f+g*u)-T-s1*g1+h'*Q*h+2*h'*S*u+u'*R*u;
    %cond2 = V-beta*(e1^2+e2^2+e3^2+e4^2+e5^2)^nv;
    %cond3 = T-beta*(e1^2+e2^2+e3^2+e4^2+e5^2)^nt;
    cond4 = R0-R;
    cond5 = S*inv(R0)*S0'+S0*inv(R0)*S'-2*S0*inv(R0)*S0'+Q0-Q;
    prog = sosineq(prog, cond1);
    %prog = sosineq(prog, cond2);
    %prog = sosineq(prog, cond3);
    %prog = sosineq(prog, R);
    prog = sosineq(prog, cond4);
    prog = sosineq(prog, cond5);


    sol = sossolve(prog);

    Q = double(sosgetsol(sol,Q));
    R = double(sosgetsol(sol,R));
    S = double(sosgetsol(sol,S));

    Q0 = Q;
    S0 = S;
    R0 = R;

    d = S*inv(R)*S'-Q;
    min(eig(d))
    i=i+1

    if min(eig(d)) >= 0 
        K = -inv(R)*S';
        V = sosgetsol(sol,V);
        break
    end
end
end
