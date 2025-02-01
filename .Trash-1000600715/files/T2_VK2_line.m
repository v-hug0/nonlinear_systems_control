clc
clear all

pvar x1 x2 x3;
x = [x1; x2; x3];

f = [x2;
    9.8-28.64*(x3+0.585)^2;
    -322.58*x1*(x3+0.585)];

g = [0;0;0.13];

% V inicial 
V = x'*x;


i = 0;
i_max = 10;
rho = 0.1;
rho_i = rho;
%alpha = 1e-6*(x'*x)^3;
while i<i_max
    %%%%%%
    prog = sosprogram(x);
    %.....controller
    [prog,u] = sospolyvar(prog,[monomials(x,0:2)],'wscoeff');
    [prog,alfa] = sospolyvar(prog,[monomials(x,2:6)],'wscoeff');
    [prog,l2] = sospolyvar(prog,[monomials(x,2:6)],'wscoeff');
    %.....SOS constraints
    gradV = [diff(V,x1) diff(V,x2) diff(V,x3)];
    dV = -gradV*(f+g*u)-alfa*(rho-V)-l2; % dV after S-procedure
    prog = sosineq(prog,dV);

    sol = sossolve(prog);
    u = sosgetsol(sol,u);
    alfa = sosgetsol(sol,alfa);
    l2 = sosgetsol(sol,l2);
 
    %%%%%%%
    dpvar rho;
    prog = sosprogram(x);
    prog = sosdecvar(prog,rho);
    [prog,V] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    [prog,l1] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    [prog,l2] = sospolyvar(prog,[monomials(x,2:6)],'wscoeff');
    %[prog,l2] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
    %.....SOS constraints
    gradV = [diff(V,x1) diff(V,x2) diff(V,x3)];
    dV = -gradV*(f+g*u)-alfa*(rho-V)-l2;
    prog = sosineq(prog,dV); 
    prog = sosineq(prog,V-l1);
    prog = sosineq(prog,rho-rho_i);
    %prog = sosineq(prog,-rho+1);
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    rho = sosgetsol(sol,rho)
    rho_i=rho;
    %l2 = sosgetsol(sol,l2);
    i = i+1
        % Feas perto de 1
    if abs(1-sol.solinfo.info.feasratio)<=0.2    
        break
    end
end

%%

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % RESPOSTA EM MALHA FECHADA
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [1;-0.7;-0.4];
[t, x] = ode45(@closedloop, [0 15], x0);
response_CL = plot(t,x(:,1),t,x(:,2),t,x(:,3));
title('Malha Fechada')
legend('x1','x2','x3','Location','best')

function dx = closedloop(t,x)
    % para l2=4, e u=0:2
    k = -7.7919*x(1)^2 + 2.3127*x(1)*x(2) + 163.2668*x(1)*(x(3)-0.585) + 11.2868*x(2)^2 + 15.9523*x(2)*(x(3)-0.585) - 9.5115*(x(3)-0.585)^2 + 109.6672*x(1) - 32.1088*x(2) - 449.3299*(x(3)-0.585) - 0.0024542
    % para l2=6, e u=0:2
    %k =   1350.3279*x(1)^2 + 561.3422*x(1)*x(2) + 261.7106*x(1)*(x(3)-0.585) - 90.2592*x(2)^2 - 260.6939*x(2)*(x(3)-0.585) + 197.0053*(x(3)-0.585)^2 + 3386.9119*x(1) + 27.8334*x(2) - 914.8812*(x(3)-0.585) + 0.0038421
    dx = [x(2);
          9.8-28.64*(x(3)^2;
         -322.58*x(1)*(x(3)-0.585) + 0.13*k];
end
