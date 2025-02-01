%Oscilador de Van der Pol

% Malha aberta, u=0
syms x1 x2
a=0.8;
[x1, x2] = meshgrid(-a*5:1e-1:a*5, -a*5:1e-1:a*5);
dx1 = x2;
dx2 = -x1 + (1 - x1.^2).*x2;



figure
streamslice(x1,x2,dx1,dx2,4)



% Controle não linear u(x) 
syms x1 x2
a=2;
[x1, x2] = meshgrid(-a*5:1e-1:a*5, -a*5:1e-1:a*5);
u = -1.2681*x1.^3+0.19992*x1.^2.*x2+0.030579*x1.*x2.^2-1.8445*x2.^3+4.4664*x1.^2- .2225*x1.*x2-0.64606*x2.^2-2.1293*x1-49.3754*x2;
dx1 = x2;
dx2 = -x1 + (1-x1.^2).*x2 + u;
figure
streamslice(x1,x2,dx1,dx2,4)