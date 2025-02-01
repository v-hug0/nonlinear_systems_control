clc
clear all


pvar e1 e2 e3 e4 e5;
e=[e1;e2;e3;e4;e5];
Ur = [0.2 0.2];
 
f = [Ur(1)*(e5+1);Ur(1)*e4;Ur(2);Ur(2)*(e5+1);-Ur(2)*e4];
g = [-1 e2; 0 e1; 0 -1; 0 -e5-1; 0 e4];

% V inicial 

%V = beta*e'*e;
P = eye(length(e));
%V = e'*P*e;
V = e'*P*e;
rho =  1;
ubar = [0.2;0.2];

nv = 2;
nt = 2;

lambda = 20;
%%
i = 0;
i_max = 5;
while i<i_max
    %%%%%% solving with (V,lambda fixed)
    prog = sosprogram(e);
    %.....controller
    [prog,u1] = sospolyvar(prog,[monomials(e,1:2)],'wscoeff');
    [prog,u2] = sospolyvar(prog,[monomials(e,1:2)],'wscoeff');
    u=[u1;2];
    %.....S-procedure multipliers
    [prog,s1] = sospolyvar(prog, [monomials(e,1:4)],'wscoeff');
    [prog,alfa] = sospolyvar(prog,[monomials(e,2:6)],'wscoeff');
    %.....ND maker for dV term
    [prog,T] = sospolyvar(prog,[monomials(e,1:4)],'wscoeff');
    %.....SOS constraints
    g1 = e4^2+e5^2-1;  % estados amarrados
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    %...
    dV = -gradV*(f+g*(u+ubar))-s1*g1-lambda*T-alfa*(rho-V); % dV after S-procedure
    %V_PD = V-1e-8*(e1^2+e2^2+e3^2+e4^2+e5^2)^nv; % makes V (PD)
    %T_PD = T-1e-8*(e1^2+e2^2+e3^2+e4^2+e5^2)^nt; % makes T (PD)
    %...
    %prog = sosineq(prog,V_PD);
    %prog = sosineq(prog,T_PD);
    prog = sosineq(prog,dV);
    % first solution (T,k,alfa)
    sol = sossolve(prog);
    u = sosgetsol(sol,u);
    T = sosgetsol(sol,T);
    alfa = sosgetsol(sol,alfa);
    s1 = sosgetsol(sol,s1);
    %%%%%% solving with (T,k,alfa fixed)
    %dpvar lambda;
    prog = sosprogram(e);
    %prog = sosdecvar(prog,lambda);
    %.....LF
    [prog,V] = sospolyvar(prog,[monomials(e,1:4)],'wscoeff');
    %.....S-procedure
    [prog,s1] = sospolyvar(prog, [monomials(e,1:4)],'wscoeff');
    %.....SOS constraints
    g1 = e4^2+e5^2-1;
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    %...
    dV = -gradV*(f+g*(u+ubar))-s1*g1-lambda*T-alfa*(rho-V);
    %V_PD = V-1e-8*(e1^2+e2^2+e3^2+e4^2+e5^2)^nv;
    %T_PD = T-1e-8*(e1^2+e2^2+e3^2+e4^2+e5^2)^nt;
    %...
    %prog = sosineq(prog,V_PD);
    %prog = sosineq(prog,T_PD);
    prog = sosineq(prog,dV);
    %prog = sosineq(prog,lambda+1e-6);
    %.....optimization 
    %prog = sossetobj(prog,lambda);
    % second solution (V)
    s1 = sosgetsol(sol,s1);
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    %lambda = sosgetsol(sol,lambda)
    i = i+1
end