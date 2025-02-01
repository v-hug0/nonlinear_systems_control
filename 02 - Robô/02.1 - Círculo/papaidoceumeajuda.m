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
beta = 1e-8;

nv = 2;
nt = 2;
ubar = [0.2;0.2];
%%
i = 0;
i_max = 5;
while i<i_max
    %%%%%% solving with (V fixed)
    prog = sosprogram(e);
    %.....controller
    [prog,u1] = sospolyvar(prog,[monomials(e,1:2)],'wscoeff');
    [prog,u2] = sospolyvar(prog,[monomials(e,1:2)],'wscoeff');
    u=[u1;u2];
    %.....S-procedure multipliers
    [prog,s1] = sospolyvar(prog, [monomials(e,1:4)],'wscoeff');
    [prog,alfa] = sospolyvar(prog,[monomials(e,2:6)],'wscoeff');
    [prog,alfaV] = sospolyvar(prog,[monomials(e,2:4)],'wscoeff');
    [prog,alfaT] = sospolyvar(prog,[monomials(e,2:4)],'wscoeff');
    %.....ND maker for dV term
    [prog,T] = sospolyvar(prog,[monomials(e,1:6)],'wscoeff');
    %.....SOS constraints
    g1 = e4^2+e5^2-1;  % estados amarrados
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    %...
    dV = -gradV*(f+g*(u+ubar))-s1*g1-T-alfa*(rho-V); % dV after S-procedure
    T_PD = T-beta*(e1^2+e2^2+e3^2+e4^2+e5^2)^nt-alfaT*(1-V); % makes T (PD)
    %...
    prog = sosineq(prog,T_PD);
    prog = sosineq(prog,dV);
    % first solution (T,u,alfa)
    sol = sossolve(prog);
    u = sosgetsol(sol,u);
    T = sosgetsol(sol,T);
    alfa = sosgetsol(sol,alfa);
    alfaV = sosgetsol(sol,alfaV);
    alfaT = sosgetsol(sol,alfaT);
    s1 = sosgetsol(sol,s1);
    %%%%%% solving with (u,alfa fixed)
    prog = sosprogram(e);
    %.....LF
    [prog,V] = sospolyvar(prog,[monomials(e,1:4)],'wscoeff');
    %.....ND maker for dV term
    [prog,T] = sospolyvar(prog,[monomials(e,2:6)],'wscoeff');
    %.....S-procedure
    [prog,s1] = sospolyvar(prog, [monomials(e,1:4)],'wscoeff');
    %.....SOS constraints
    g1 = e4^2+e5^2-1;
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    %...
    dV = -gradV*(f+g*(u+ubar))-s1*g1-T-alfa*(rho-V);
    V_PD = V-beta*(e1^2+e2^2+e3^2+e4^2+e5^2)^nv-alfaV*(1-V);
    T_PD = T-beta*(e1^2+e2^2+e3^2+e4^2+e5^2)^nt-alfaT*(1-V);
    %...
    prog = sosineq(prog,V_PD);
    prog = sosineq(prog,T_PD);
    prog = sosineq(prog,dV);
    % second solution (V)
    s1 = sosgetsol(sol,s1);
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    i = i+1
end
