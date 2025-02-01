clc
clear all


pvar e1 e2 e3 e4 e5;
e=[e1;e2;e3;e4;e5];
Ur = [0.2 0.2];
 
f = [Ur(1)*(e5+1);Ur(1)*e4;Ur(2);Ur(2)*(e5+1);-Ur(2)*e4];
g = [-1 e2; 0 e1; 0 -1; 0 -e5-1; 0 e4];

% V inicial 
%beta = 1;
%V = beta*e'*e;
P = eye(length(e));
%V = e'*P*e;
V = e'*P*e;
rho =  1;
lambda = 50;

%%
i = 0;
i_max = 5;
while i<i_max
    % solving with SOSTOOLS (V and ro fixed)
    prog = sosprogram(e);
    [prog,kx1] = sospolyvar(prog,[monomials(e,1:2)],'wscoeff');
    [prog,kx2] = sospolyvar(prog,[monomials(e,1:2)],'wscoeff');
    [prog,T] = sospolyvar(prog,[monomials(e,2:6)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e,1:4)],'wscoeff');
    [prog,alfa] = sospolyvar(prog,[monomials(e,2:6)],'wscoeff');
    [prog,L] = sospolyvar(prog,[monomials(e,2:4)],'wscoeff');
    kx=[kx1;kx2];
    % SOS constraints
    g1 = e4^2+e5^2-1;
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*kx)-s1*g1-T-alfa*(rho-V)+lambda*L;
    prog = sosineq(prog,T-1e-6*(e1^2+e2^2+e3^2+e4^2+e5^2));
    prog = sosineq(prog,dV);
    % first solution (k and alfa)
    sol = sossolve(prog);
    kx = sosgetsol(sol,kx);
    alfa = sosgetsol(sol,alfa);
    s1 = sosgetsol(sol,s1);
    L = sosgetsol(sol,L);
    T = sosgetsol(sol,T);
    % solving with SOSTOOLS (k,alfa and L fixed)
    dpvar lambda;
    prog = sosprogram(e);
    prog = sosdecvar(prog,lambda);
    [prog,V] = sospolyvar(prog,[monomials(e,2:4)],'wscoeff');
    [prog,T] = sospolyvar(prog,[monomials(e,2:6)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e,1:4)],'wscoeff');
    % SOS constraints
    prog = sosineq(prog,V);
    g1 = e4^2+e5^2-1;
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*kx)-s1*g1-T-alfa*(rho-V)+lambda*L;
    prog = sosineq(prog,V-1e-6*(e1^2+e2^2+e3^2+e4^2+e5^2));
    prog = sosineq(prog,T-1e-6*(e1^2+e2^2+e3^2+e4^2+e5^2));
    prog = sosineq(prog,dV);
    prog = sossetobj(prog,lambda);
    % second solution (V and lambda)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    s1 = sosgetsol(sol,s1);
    lambda = sosgetsol(sol,lambda);
    T = sosgetsol(sol,T);
    if lambda<0
        break
    else
        i = i+1
    end
end