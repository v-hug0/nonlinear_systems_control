clear all
close all
clc

%%
% Inicializações - CIRCULO
r = 1;
w = 2;
syms t
x = r*cos(w*t);
y = r*sin(w*t);
dx = diff(x,t);
dy = diff(y,t);
ddx = diff(dx,t);
ddy = diff(dy,t);

global vr wr dvr dwr

vr = sqrt(dx^2+dy^2);
wr = (dx*ddy-dy*ddx)/(dx^2+dy^2);
dvr = diff(vr,t);
dwr = diff(wr,t);
vr = matlabFunction(vr);
wr = matlabFunction(wr);
dvr = matlabFunction(dvr);
dwr = matlabFunction(dwr);

%%

% Inicializações - LEMNISCATA
l = 2.5; % 2500 mm
w = pi/10;
syms t
x = l*sin(2*w*t);
y = l*cos(w*t);

dx = diff(x,t);
dy = diff(y,t);

ddx = diff(dx,t);
ddy = diff(dy,t);

global vr wr dvr dwr

vr = sqrt(dx^2+dy^2);
wr = (dx*ddy-dy*ddx)/(dx^2+dy^2);
dvr = diff(vr,t);
dwr = diff(wr,t);
vr = matlabFunction(vr);
wr = matlabFunction(wr);
dvr = matlabFunction(dvr);
dwr = matlabFunction(dwr);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           TRAJECTORY GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inicializations for reference trajectory
dt = 0.1;          % sampling time
Tsim = 10;   %90      % tempo de simulacao
h = dt;           % integrating time
xr = [0;-1;0];
Xr = xr;         %Euler


% Inicializations for real trajectory

%xo = [-1;-2.5;0];
xo = [-0.3;-0.8;0];
x = xo;
X = x;
v = [];
w = [];
syms E1 E2 E4 E5
E = [E1;E2;E4;E5-1];
Z = monomials(E,1:4);



for k = 1:round(Tsim/h)  % 1 a tempo de simulacao/ tempo de integracao

    %%%%% SOSPROGRAM FOR EACH POINT
    pvar e1 e2 e4 e5
    e_n=[e1;e2;e4;e5];
    ubar = [vr(k*dt),wr(k*dt)];
    f = [vr(k*dt)*e5+e2*ubar(2); vr(k*dt)*e4-e1*ubar(2); wr(k*dt)*e5; 0];
    g = [-1 e2; 0 -e1; 0 -e5-1; 0 e4];
    
    % Initial V to be fixed in the first program of the iteration 1.
    Vo = e_n'*eye(length(e_n))*e_n;
    V = Vo;
    
    % Maximum iteration.
    i = 0;
    i_max = 100;
    
    % Monomials of the controllers [v;ω] whose each coefficient (the controller
    % gains) will be found as part of the SOS Program solution.
    z = monomials(e_n,1:4);
    
    % Tolerance β between the previous and current solution
    beta = 0.5;
    
    % Inicial choice for the controller coefficients to enable the first check
    % of the stopping criterion. We do set a large value to all coefficient so the first check
    % be guaranteed to not fall within the tolerance β.
    uo = ones(length(monomials(e_n,1:4)),size(g,2))*100;
    
    % Factibility loop until solution convergence.
    while i<i_max
        %%%%% SOS-PROGRAM 1: SET OF CONTROLLER GAINS SOLUTION FROM THE PREVIOUS FIXED V %%%%%
        % solving with SOSTOOLS (V fixed).
        % Defining the first set of decision variables (u,l2,s1).
        prog = sosprogram(e_n);
        [prog,K1] = sospolymatrixvar(prog,monomials(e_n,0),[1,length(z)]);
        [prog,K2] = sospolymatrixvar(prog,monomials(e_n,0),[1,length(z)]);
        u=[K1*z;K2*z];
        [prog,l2] = sossosvar(prog,[monomials(e_n,2:4)],'wscoeff');
        [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)],'wscoeff');
        % Defining the SOS constraints
        g1 = e4^2+e5^2-1; % constraint from eq.(13)
        gradV = [diff(V,e1);diff(V,e2);diff(V,e4);diff(V,e5)];
        dV = -gradV'*(f+g*u)-s1*g1-l2;
        prog = sosineq(prog,dV);
        % First solution: (u,l2,s1)
        sol = sossolve(prog);
        u = sosgetsol(sol,u);
        K1 = double(sosgetsol(sol,K1));
        K2 = double(sosgetsol(sol,K2));
        l2 = sosgetsol(sol,l2);
        s1 = sosgetsol(sol,s1);
        %%%%% SOS-PROGRAM 2: LYAPUNOV FUNCTION SOLUTION FROM THE PREVIOUS CONTROLLER GAINS %%%%%
        % solving with SOSTOOLS (u fixed)
        % Defining the second set of decision variables (V,,l1,l2,s1).
        prog = sosprogram(e_n);
        [prog,V] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
        [prog,l1] = sossosvar(prog,[monomials(e_n,2:4)],'wscoeff');
        [prog,l2] = sossosvar(prog,[monomials(e_n,2:4)],'wscoeff');
        [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)], 'wscoeff');
        %[prog,s2] = sospolyvar(prog, [monomials(e_n,2:4)], 'wscoeff');
        % Defining the SOS constraints
        g1 = e4^2+e5^2-1; % constraint from eq.(13)
        gradV = [diff(V,e1);diff(V,e2);diff(V,e4);diff(V,e5)];
        dV = -gradV'*(f+g*u)-s1*g1-l2;
        prog = sosineq(prog,dV);
        prog = sosineq(prog,V-l1);
        % Second solution: (V,l1,l2,s1)
        sol = sossolve(prog);
        V = sosgetsol(sol,V);
        l1 = sosgetsol(sol,l1);
        l2 = sosgetsol(sol,l2);
        s1 = sosgetsol(sol,s1);
        %s2 = sosgetsol(sol,s2);
        % Stopping criterion
        delta = max(abs(full(u.coefficient)-uo))
        if delta<beta
            break
        end
        % Stores the current controller coefficients to be compared in the next
        % loop.
        uo = full(u.coefficient);
        i = i+1
    end

    K11(k,:) = K1;
    K22(k,:) = K2;

    %%%%%%% 
    Xr = [Xr xr+h*([vr(k*dt)*cos(xr(3));vr(k*dt)*sin(xr(3));wr(k*dt)])];  %Euler
    xr = Xr(:,k+1);
    plot(Xr(1,:),Xr(2,:),'b.'); drawnow
    hold on
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(2);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);
    v(k) = K11(k,:)*subs(Z)+ubar(1);
    w(k) = K22(k,:)*subs(Z)+ubar(2);
    X = [X x+dt*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
    plot(X(1,:),X(2,:),'+'); drawnow
end

%%
figure
%hold on
plot(Xr(1,:),Xr(2,:),'');
hold on
plot(X(1,:),X(2,:),'');


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          ROBOT REFERENCE TRACKING SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [-1;-2.5;0.9];
x = xo;
X = x;
v = [];
w = [];

syms E1 E2 E4 E5
E = [E1;E2;E4;E5-1];
Z = monomials(E,1:4);

for k = 1:(round(Tsim/h))    % 1 a tempo de simulacao/tempo de integracao  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(2);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);

    %v(k) = K1*subs(Z)+ubar(1);
    %w(k) = K2*subs(Z)+ubar(2);
    v(k) = -0.075303*E1^4 - 0.030166*E1^3*E2 - 0.087852*E1^3*E4 - 0.27804*E1^3*(E5-1) + 0.028573*E1^2*E2^2 + 0.10772*E1^2*E2*E4 + 0.16195*E1^2*E2*(E5-1) - 0.11447*E1^2*E4^2 - 0.21984*E1^2*E4*(E5-1) - 0.12463*E1^2*(E5-1)^2 + 0.0043471*E1*E2^3 - 0.046854*E1*E2^2*E4 - 0.15103*E1*E2^2*(E5-1) - 0.0070795*E1*E2*E4^2 - 0.0047084*E1*E2*E4*(E5-1) + 0.062875*E1*E2*(E5-1)^2 - 0.012684*E1*E4^3 - 0.1264*E1*E4^2*(E5-1) - 0.026643*E1*E4*(E5-1)^2 - 0.24807*E1*(E5-1)^3 + 0.07716*E2^4 + 0.11832*E2^3*E4 + 0.14263*E2^3*(E5-1) + 0.056355*E2^2*E4^2 - 0.2546*E2^2*E4*(E5-1) - 0.03893*E2^2*(E5-1)^2 + 0.2447*E2*E4^3 + 0.45614*E2*E4^2*(E5-1) + 0.52292*E2*E4*(E5-1)^2 + 0.29669*E2*(E5-1)^3 + 0.10717*E4^4 + 0.27471*E4^3*(E5-1) + 0.33554*E4^2*(E5-1)^2 + 0.18185*E4*(E5-1)^3 + 0.020531*(E5-1)^4 + 0.38561*E1^3 - 1.1464*E1^2*E2 - 1.9633*E1^2*E4 - 3.6263*E1^2*(E5-1) + 0.56878*E1*E2^2 - 0.64452*E1*E2*E4 - 0.94175*E1*E2*(E5-1) + 0.069601*E1*E4^2 - 1.4704*E1*E4*(E5-1) - 0.2084*E1*(E5-1)^2 - 0.99687*E2^3 - 1.6795*E2^2*E4 - 3.05*E2^2*(E5-1) - 1.602*E2*E4^2 - 2.4894*E2*E4*(E5-1) - 2.5062*E2*(E5-1)^2 - 1.7428*E4^3 - 3.6391*E4^2*(E5-1) - 2.8992*E4*(E5-1)^2 - 3.7739*(E5-1)^3 - 0.076835*E1^2 + 0.07431*E1*E2 - 0.48516*E1*E4 + 0.050438*E1*(E5-1) + 0.9536*E2^2 + 0.33252*E2*E4 + 4.5428*E2*(E5-1) - 3.6971*E4^2 - 8.8372*E4*(E5-1) + 0.93152*(E5-1)^2 + 2.3975*E1 - 1.6345*E2 - 2.4543*E4 + 0.63279*(E5-1) + ubar(1);
    w(k) = -0.16276*E1^4 - 0.034436*E1^3*E2 - 0.24106*E1^3*E4 - 0.82512*E1^3*(E5-1) + 0.067252*E1^2*E2^2 + 0.42547*E1^2*E2*E4 + 0.71967*E1^2*E2*(E5-1) - 0.14095*E1^2*E4^2 - 0.25793*E1^2*E4*(E5-1) - 0.65219*E1^2*(E5-1)^2 - 0.075421*E1*E2^3 - 0.33152*E1*E2^2*E4 - 0.81321*E1*E2^2*(E5-1) - 0.12481*E1*E2*E4^2 - 0.036941*E1*E2*E4*(E5-1) + 0.16753*E1*E2*(E5-1)^2 - 0.21964*E1*E4^3 - 0.64916*E1*E4^2*(E5-1) - 0.10142*E1*E4*(E5-1)^2 - 0.84357*E1*(E5-1)^3 + 0.22029*E2^4 + 0.4452*E2^3*E4 + 0.6*E2^3*(E5-1) + 0.26459*E2^2*E4^2 - 0.34241*E2^2*E4*(E5-1) - 0.25408*E2^2*(E5-1)^2 + 0.43711*E2*E4^3 + 0.53234*E2*E4^2*(E5-1) + 0.374*E2*E4*(E5-1)^2 + 0.68238*E2*(E5-1)^3 + 0.052412*E4^4 - 0.16557*E4^3*(E5-1) - 0.25174*E4^2*(E5-1)^2 - 0.20562*E4*(E5-1)^3 - 0.43588*(E5-1)^4 - 0.15071*E1^3 - 0.071845*E1^2*E2 - 0.23153*E1^2*E4 - 0.8718*E1^2*(E5-1) - 0.083079*E1*E2^2 - 0.005279*E1*E2*E4 + 0.098609*E1*E2*(E5-1) - 0.21204*E1*E4^2 - 0.16585*E1*E4*(E5-1) - 0.10272*E1*(E5-1)^2 - 0.057006*E2^3 - 0.09572*E2^2*E4 - 0.59361*E2^2*(E5-1) - 0.11159*E2*E4^2 - 0.027012*E2*E4*(E5-1) - 0.038652*E2*(E5-1)^2 - 0.23668*E4^3 - 0.75296*E4^2*(E5-1) + 0.37186*E4*(E5-1)^2 - 0.52779*(E5-1)^3 - 0.12877*E1^2 + 0.039287*E1*E2 - 0.59691*E1*E4 - 0.93698*E1*(E5-1) - 0.04627*E2^2 - 0.48294*E2*E4 - 0.40361*E2*(E5-1) - 1.2636*E4^2 - 4.1331*E4*(E5-1) - 1.969*(E5-1)^2 - 0.18379*E1 + 1.2361*E2 + 2.0829*E4 + 4.7798*(E5-1) + ubar(2);

    X = [X x+dt*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
    
end

Error = Xr-X;
Error4 = sin(Error(3,:));
Error5 = cos(Error(3,:));

plot(X(1,:),X(2,:),'r')
legend('Reference Robot','Real Robot')
xlabel('x (m)')
ylabel('y (m)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem de v e w 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
time = dt*(0:Tsim/h-1);
subplot(2,1,1)
plot(time,v,'b')
xlabel('Time (s)')
ylabel('U1 (m/s)')
title('v')

hold on
subplot(2,1,2)
plot (time,w,'r')
xlabel('Time (s)')
ylabel('U2 (rad/s)')
title('ω')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem do erro x y e theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(time,Error(1,1:end-1),time,Error(2,1:end-1),time,Error(3,1:end-1),time,Error4(1:end-1),time,Error5(1:end-1));
legend('e1','e2','e3','e4','e5')
