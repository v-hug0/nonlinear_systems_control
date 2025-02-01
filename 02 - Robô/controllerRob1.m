
clear all
close all
clc

pvar e1 e2 e3 e4 e5
e_n=[e1;e2;e3;e4;e5];
Ur = [0.1 0];
f = [Ur(1)*(e5+1); Ur(1)*e4; Ur(2); Ur(2)*(e5+1); -Ur(2)*e4];
g = [-1 e2; 0 -e1; 0 -1; 0 -e5-1; 0 e4];

n = size(f,2);
m = size(g,2);

% V inicial 
V = e_n'*e_n;
i = 0;
i_max = 20;
while i<i_max
    % solving with SOSTOOLS (V fixed)
    prog = sosprogram(e_n);
    [prog,kx1] = sospolyvar(prog,[monomials([e1;e2],1:2)],'wscoeff');
    [prog,kx2] = sospolyvar(prog,[monomials([e2;e3;e4;e5],1:2)],'wscoeff');
    kx=[kx1;kx2];
    [prog,T] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,1:4)],'wscoeff');
    g1 = e4^2+e5^2-1;
    % SOS constraints
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*kx)-T-s1*g1;
    prog = sosineq(prog,dV);
    % first solution (k and lambda)
    sol = sossolve(prog);
    kx = sosgetsol(sol,kx);
    % solving with SOSTOOLS (k fixed)
    prog = sosprogram(e_n);
    [prog,V] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,1:4)], 'wscoeff');
    [prog,T] = sospolyvar(prog,[monomials(e_n,1:4)],'wscoeff');
    g1 = e4^2+e5^2-1;
    % SOS constraints
    prog = sosineq(prog,V);
    gradV = [diff(V,e1) diff(V,e2) diff(V,e3) diff(V,e4) diff(V,e5)];
    dV = -gradV*(f+g*kx)-T-s1*g1;
    prog = sosineq(prog,dV);
    % first solution (V)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    i = i+1
end



