clear
close all
clc
syms x1 x2 x3;
var = [x1;x2;x3];
e = 1;

f = [-x1^3-x1*x3^2; -x2-x1^2*x2; -x3-3*x3/(x3^2+1)+3*x1^2*x3];

prog = sosprogram(var);
[prog, V] = sospolyvar(prog,[x1^2;x2^2;x3^2],'wscoeff');
prog = sosineq(prog,V-e*(x1^2 + x2^2 + x3^2));
prog = sosineq(prog,-(x3^2+1)*(diff(V,x1)*f(1)+diff(V,x2)*f(2)+diff(V,x3)*f(3)));
sol = sossolve(prog);
SOLV = sosgetsol(sol,V)

a = 5;
[x1,x2,x3] = meshgrid(-a:1:a,-a:1:a,-a:1:a);
[x1,x3] = meshgrid(-a:1:a,-a:1:a);
dx1 = -x1.^3-x1.*x3.^2;
% dx2 = -x2-x1.^2.*x2;
dx3 = -x3-3.*x3./(x3.^2+1)+3.*x1.^2.*x3;
streamslice(x1,x3,dx1,dx3)