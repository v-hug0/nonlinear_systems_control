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

pvar e1 e2 e4 e5 u1 u2
var = [e1;e2;e4;e5;u1;u2];
e_n = [e1;e2;e4;e5];
u = [u1;u2];
Ur = [0.2;0.2];
ubar = Ur;
f = [Ur(1)*(e5+1);
     Ur(1)*e4;
     Ur(2)*(e5+1);
    -Ur(2)*e4];
g = [-1 e2;
      0 -e1;
      0 -e5-1;
      0 e4];

% h1 = monomials(e_n,1);
% h1 = ones(1,length(h1))*h1;
% h2 = monomials(e_n,2);
% h2 = ones(1,length(h2))*h2;
% h3 = monomials(e_n,3);
% h3 = ones(1,length(h3))*h3;
% h4 = monomials(e_n,4);
% h4 = ones(1,length(h4))*h4;
% h = [h1;h2;h3;h4];
h = monomials(e_n,1:4);

n = size(f,2);
m = size(g,2);
p = size(h,1); %numero de linhas do h
j = size(h,2);

i = 0;
imax = 100; 
%beta = 0.01; 
z = 12;

prog=sosprogram(var);
[prog,V0] = sospolyvar(prog,monomials(e_n,2:4),'wscoeff');
[prog,l1] = sospolyvar(prog,monomials(e_n,2:4),'wscoeff');
[prog,l2] = sospolyvar(prog,monomials(e_n,2:4),'wscoeff');
[prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)],'wscoeff');
[prog,q]=sospolymatrixvar(prog,monomials(var,0),[z,z],'symmetric');
[prog,I]=sospolymatrixvar(prog,monomials(var,0),[p-z,p-z],'symmetric');
Q = [q zeros(z,p-z);
     zeros(p-z,z) I*eye(p-z)];
[prog,S]=sospolymatrixvar(prog,monomials(var,0),[p,m]);
[prog,R]=sospolymatrixvar(prog,monomials(var,0),[m,m]);
%R = [R(1) 0; 0 R(2)];
%[prog,R]=sospolymatrixvar(prog,monomials(var,0),[m,m],'symmetric');
g1 = e4^2+e5^2-1;

gradV0 = [diff(V0,e1) diff(V0,e2) diff(V0,e4) diff(V0,e5)];
cond1 = -gradV0*(f+g*(u+ubar))+h'*Q*h+2*h'*S*(u+ubar)+(u+ubar)'*R*(u+ubar)-s1*g1-l2;
% cond1 = -[gradV0*f-h'*Q*h+l2 1/2*gradV0*g-h'*S;
%           (1/2*gradV0*g-h'*S)' -R];
cond2 = V0-l1;
prog = sosineq(prog, cond1);
prog = sosineq(prog, cond2);
prog = sosineq(prog, R);

sol = sossolve(prog);

Q0 = double(sosgetsol(sol,Q));
R0 = double(sosgetsol(sol,R));
S0 = double(sosgetsol(sol,S));


d = S0*inv(R0)*S0'-Q0;

if d >=0
    K = -inv(R0)*S0';
    V0 = sosgetsol(sol,V0);
else 
    while i<imax
    prog=sosprogram(var);
    [prog,V] = sospolyvar(prog,monomials(e_n,2:4),'wscoeff');
    [prog,l1] = sospolyvar(prog,monomials(e_n,2:4),'wscoeff');
    [prog,l2] = sospolyvar(prog,monomials(e_n,2:4),'wscoeff');
    [prog,s1] = sospolyvar(prog,monomials(e_n,2:4),'wscoeff');
    [prog,I]=sospolymatrixvar(prog,monomials(var,0),[p-z,p-z],'symmetric');
    Q = [q zeros(z,p-z);
         zeros(p-z,z) I*eye(p-z)];;
    [prog,S]=sospolymatrixvar(prog,monomials(var,0),[p,m]); 
    [prog,R]=sospolymatrixvar(prog,monomials(var,0),[m,m],'symmetric'); 
    g1 = e4^2+e5^2-1;
    
    gradV = [diff(V,e1) diff(V,e2) diff(V,e4) diff(V,e5)];
    cond1 = -gradV*(f+g*(u+ubar))+h'*Q*h+2*h'*S*(u+ubar)+(u+ubar)'*R*(u+ubar)-s1*g1-l2;
    % cond1 = -[gradV*f-h'*Q*h+l2 1/2*gradV*g-h'*S;
    %       (1/2*gradV*g-h'*S)' -R];
    cond2 = V-l1;
    cond5 = R0-R;
    cond6 = S*inv(R0)*S0'+S0*inv(R0)*S'-2*S0*inv(R0)*S0'+Q0-Q;
    prog = sosineq(prog, cond1);
    prog = sosineq(prog, cond2);
    prog = sosineq(prog, R);
    prog = sosineq(prog, cond5);
    prog = sosineq(prog, cond6);

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
