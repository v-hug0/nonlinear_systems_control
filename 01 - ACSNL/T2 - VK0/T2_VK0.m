clc
clear all

pvar x1 x2 x3;
x = [x1; x2; x3];


M = 653;
N = 324;
A = 0.0235;
R = 0.5;
G = 9.8;
mi = 4*pi*1e-7;

%xo 


f = [x2;
     G-A/(mi*M)*x3^2;
     -(2*R)/(mi*A*N^2)*x1*x3];

g = [0;0;1/(N*A)];

% V inicial 
V = x'*x;

i = 0;
i_max = 20;
while i<i_max
    %%%%%%
    prog = sosprogram(x);
    %.....controller
    [prog,u] = sospolyvar(prog,[monomials(x,0:2)],'wscoeff');
    [prog,alfa] = sospolyvar(prog,[monomials(x,2:6)],'wscoeff');
    [prog,l2] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    %.....SOS constraints
    gradV = [diff(V,x1) diff(V,x2) diff(V,x3)];
    dV = gradV*(f+g*u)+alfa*(1-V)+l2; % dV after S-procedure
    prog = sosineq(prog,-dV);

    sol = sossolve(prog);
    u = sosgetsol(sol,u);
    alfa = sosgetsol(sol,alfa);
    l2 = sosgetsol(sol,l2);

    %%%%%%%
    prog = sosprogram(x);
    [prog,V] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    [prog,l1] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    [prog,l2] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    %.....SOS constraints
    gradV = [diff(V,x1) diff(V,x2) diff(V,x3)];
    dV = gradV*(f+g*u)+alfa*(1-V)+l2;
    prog = sosineq(prog,-dV); 
    prog = sosineq(prog,V-l1);

    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    l2 = sosgetsol(sol,l2);
    i = i+1
end

%%

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % RESPOSTA EM MALHA FECHADA
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [0.008;-0.003;-0.4];
[t, x] = ode45(@closedloop, [0 40], x0);
response_CL = plot(t,x(:,1),t,x(:,2),t,x(:,3));
title('Malha Fechada')
legend('x1','x2','x3','Location','best')

function dx = closedloop(t,x)
    % para l2=4, e u=0:2
    k =   1304.2899*x(1)^2 + 150.8358*x(1)*x(2) + 946.5424*x(1)*(x(3)-0.585) - 376.8834*x(2)^2 + 117.3276*x(2)*(x(3)-0.585) - 527.1092*(x(3)-0.585)^2 + 1302.9994*x(1) + 88.0924*x(2) - 1586.114*(x(3)-0.585) - 0.0032884
    % para l2=6, e u=0:2
    %k =   1350.3279*x(1)^2 + 561.3422*x(1)*x(2) + 261.7106*x(1)*(x(3)-0.585) - 90.2592*x(2)^2 - 260.6939*x(2)*(x(3)-0.585) + 197.0053*(x(3)-0.585)^2 + 3386.9119*x(1) + 27.8334*x(2) - 914.8812*(x(3)-0.585) + 0.0038421
    %k =   674.635*x(1)^2 - 183.0562*x(1)*x(2) + 450.0423*x(1)*(x(3)-0.585) - 68.1116*x(2)^2 + 179.0579*x(2)*(x(3)-0.585) + 23.7489*(x(3)-0.585)^2 + 580.9475*x(1) + 464.3459*x(2) - 952.2373*(x(3)-0.585) - 0.0084595
    %k = 399.6652*x(1)^2 + 410.8575*x(1)*x(2) + 1978.0395*x(1)*x(3) - 85.6032*x(2)^2 - 43.8547*x(2)*x(3) - 269.212*x(3)^2 + 723.2605*x(1) - 1335.7506*x(2) - 1358.9552*x(3) - 0.0077202
    dx = [x(2);
          9.8-28.64*x(3)^2;
         -322.58*x(1)*x(3) + 0.13*k];
end
