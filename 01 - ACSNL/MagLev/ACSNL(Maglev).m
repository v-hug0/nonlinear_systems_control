clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SUM-OF-SQUARES FEASIBILITY (OR OPTIMIZATION) PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initializations
pvar x1 x2 x3;
x=[x1;x2;x3];
ubar = 0;
G = 9.81;
R = 0.5;
A = 0.0235;
mi = 4*pi*1e-7;
M = 653;
N = 342;


% Polynomial recasted robot tracking model from eq.(16).
f = [x2;
     G-A/(mi*M)*x3^2;
    -2*R/(mi*A*N^2)*x1*x3];

g = [0;0;1/(N*A)];

% Initial V to be fixed in the first program of the iteration 1.
Vo = x'*eye(length(x))*x;
V = Vo;

% Maximum iteration.
i = 0;
i_max = 100;

% Monomials of the controllers [v;ω] whose each coefficient (the controller
% gains) will be found as part of the SOS Program solution.
z = monomials(x,1:4);


% Factibility loop until solution convergence.
while i<i_max
    %%%%% SOS-PROGRAM 1: SET OF CONTROLLER GAINS SOLUTION FROM THE PREVIOUS FIXED V %%%%%
    % solving with SOSTOOLS (V fixed).
    % Defining the first set of decision variables (u,l2,s1).
    prog = sosprogram(x);
    [prog,K] = sospolymatrixvar(prog,monomials(x,0),[1,length(z)]);
    u=K*z;
    [prog,l2] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(x,2:4)],'wscoeff');
    % Defining the SOS constraints
    gradV = [diff(V,x1);diff(V,x2);diff(V,x3)];
    g1 = 1-V;
    dV = -gradV'*(f+g*(u+ubar))-s1*g1-l2;
    prog = sosineq(prog,dV);
    % First solution: (u,l2,s1)
    sol = sossolve(prog);
    u = sosgetsol(sol,u);
    K = double(sosgetsol(sol,K));
    l2 = sosgetsol(sol,l2);
    s1 = sosgetsol(sol,s1);
    %%%%% SOS-PROGRAM 2: LYAPUNOV FUNCTION SOLUTION FROM THE PREVIOUS CONTROLLER GAINS %%%%%
    % solving with SOSTOOLS (u fixed)
    % Defining the second set of decision variables (V,,l1,l2,s1).
    prog = sosprogram(x);
    [prog,V] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    [prog,l1] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    [prog,l2] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    % Defining the SOS constraints
    g1 = 1-V;
    gradV = [diff(V,x1);diff(V,x2);diff(V,x3)];
    dV = -gradV'*(f+g*(u+ubar))-s1*g1-l2;
    prog = sosineq(prog,dV);
    prog = sosineq(prog,V-l1);
    % Second solution: (V,l1,l2,s1)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    l1 = sosgetsol(sol,l1);
    l2 = sosgetsol(sol,l2);

    i = i+1
end


%%

x0 = [0.5; 0.5; 0.5];

[t, x] = ode45(@closedloop, [0 20], x0);
response_CL = plot(t,x(:,1),t,x(:,2),t,x(:,3));
title('Malha Fechada')
legend('x1','x2','x3','Location','best')


syms X1 X2 X3
E = [E1;E2;X3];
Z = monomials(E,1:4);

function dx = closedloop(t,x)

    v = K*subs(Z)+ubar;

    k = - 0.60357*x(1) + 0.51918*x(2) + 2.6414*x(3) + 1.1853*x(1)^2 - 0.0039012*x(2)^2 + 0.15242*x(3)^2 - 0.2542*x(1)*x(2) - 0.58084*x(1)*x(3) + 0.15914*x(2)*x(3) + 8.0966*x(1)^3 + 0.07795*x(2)^3 + 0.25808*x(3)^3 - 1.9417*x(1)^2*x(2) - 4.2947*x(1)^2*x(3) + 0.5879*x(1)*x(2)^2 + 0.34048*x(1)*x(3)^2 + 0.077199*x(2)^2*x(3) - 0.13242*x(2)*x(3)^2 + 0.69357*x(1)*x(2)*x(3);
    dx = [  -0.878*x(1)+x(3)-x(1)^2*x(3)-0.0896*x(1)*x(3)-0.019*x(2)^2+0.473*x(1)^2+3.813*x(1)^3-0.216*k;
            x(3);
            -4.209*x(1)-0.396*x(3)-0.408*x(1)^2-2.137*x(1)^3-20.991*k];
end


