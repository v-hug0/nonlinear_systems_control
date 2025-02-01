clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SUM-OF-SQUARES FEASIBILITY (OR OPTIMIZATION) PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initializations
pvar e1 e2 e4 e5
e_n=[e1;e2;e4;e5];
Ur = [0.3;0.3];
ubar = Ur;

% Polynomial recasted robot tracking model from eq.(16).
f = [Ur(1)*(e5+1); Ur(1)*e4; Ur(2)*(e5+1); -Ur(2)*e4];
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
beta = 0.1;

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
    [prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)],'wscoeff');
    [prog,s2] = sospolyvar(prog, [monomials(e_n,2:4)], 'wscoeff');
    % Defining the SOS constraints
    g1 = e4^2+e5^2-1; % constraint from eq.(13)
    gradV = [diff(V,e1);diff(V,e2);diff(V,e4);diff(V,e5)];
    dV = -gradV'*(f+g*(u+ubar))-s1*g1-l2-s2*(1-V);
    prog = sosineq(prog,dV);
    % First solution: (u,l2,s1)
    sol = sossolve(prog);
    u = sosgetsol(sol,u);
    K1 = double(sosgetsol(sol,K1));
    K2 = double(sosgetsol(sol,K2));
    l2 = sosgetsol(sol,l2);
    s1 = sosgetsol(sol,s1);
    s2 = sosgetsol(sol,s2);
    %%%%% SOS-PROGRAM 2: LYAPUNOV FUNCTION SOLUTION FROM THE PREVIOUS CONTROLLER GAINS %%%%%
    % solving with SOSTOOLS (u fixed)
    % Defining the second set of decision variables (V,,l1,l2,s1).
    prog = sosprogram(e_n);
    [prog,V] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,l1] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)], 'wscoeff');
    % Defining the SOS constraints
    g1 = e4^2+e5^2-1; % constraint from eq.(13)
    gradV = [diff(V,e1);diff(V,e2);diff(V,e4);diff(V,e5)];
    dV = -gradV'*(f+g*(u+ubar))-s1*g1-l2-s2*(1-V);
    prog = sosineq(prog,dV);
    prog = sosineq(prog,V-l1);
    % Second solution: (V,l1,l2,s1)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    l1 = sosgetsol(sol,l1);
    l2 = sosgetsol(sol,l2);
    s1 = sosgetsol(sol,s1);

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


%delta_V = max(abs(full(V.coefficient)-Vo)) in case
%Vo = full(V.coefficient);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           REFERENCE TRAJECTORY GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Tsim = 21;
dt = 0.01;
points = Tsim/dt;

k = 1;
t = linspace(0,Tsim,points);

% CIRCULAR TRAJ.
wr = 0.7;
vr = wr;
r = 1;

thetar = pi/2+wr*t;
xr = r*cos(wr*t);
yr = r*sin(wr*t);
Xr=[xr;yr;thetar];

ubar=[0.7;0.7];



figure
plot(Xr(1,1),Xr(2,1),'o','MarkerFaceColor',	"#4DBEEE");
hold on
plot(Xr(1,:),Xr(2,:),'LineWidth',2,'Color',	"#4DBEEE");
hold on
%%

h = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          ROBOT REFERENCE TRACKING SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [-0.3;-1.3;0];
x = xo;
X = x;
v = [];
w = [];

syms E1 E2 E4 E5
E = [E1;E2;E4;E5-1];
Z = monomials(E,1:4);

for k = 1:(round(Tsim/h))    % 1 a tempo de simulacao/tempo de integracao  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(1);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);

    v(k) = -0.076708*E1^4 - 0.030329*E1^3*E2 - 0.080385*E1^3*E4 - 0.24685*E1^3*(E5-1) + 0.029066*E1^2*E2^2 + 0.10604*E1^2*E2*E4 + 0.14485*E1^2*E2*(E5-1) - 0.10566*E1^2*E4^2 - 0.1691*E1^2*E4*(E5-1) - 0.11084*E1^2*(E5-1)^2 + 0.0038438*E1*E2^3 - 0.040417*E1*E2^2*E4 - 0.13301*E1*E2^2*(E5-1) - 0.0041475*E1*E2*E4^2 + 0.0041687*E1*E2*E4*(E5-1) + 0.053506*E1*E2*(E5-1)^2 - 0.0062198*E1*E4^3 - 0.10731*E1*E4^2*(E5-1) - 0.017483*E1*E4*(E5-1)^2 - 0.21764*E1*(E5-1)^3 + 0.0775*E2^4 + 0.11462*E2^3*E4 + 0.12086*E2^3*(E5-1) + 0.065951*E2^2*E4^2 - 0.20404*E2^2*E4*(E5-1) - 0.01337*E2^2*(E5-1)^2 + 0.23784*E2*E4^3 + 0.4108*E2*E4^2*(E5-1) + 0.47556*E2*E4*(E5-1)^2 + 0.26268*E2*(E5-1)^3 + 0.11125*E4^4 + 0.27611*E4^3*(E5-1) + 0.32499*E4^2*(E5-1)^2 + 0.18826*E4*(E5-1)^3 + 0.027639*(E5-1)^4 + 0.40246*E1^3 - 1.1461*E1^2*E2 - 1.8728*E1^2*E4 - 3.1379*E1^2*(E5-1) + 0.58759*E1*E2^2 - 0.61068*E1*E2*E4 - 0.80671*E1*E2*(E5-1) + 0.095044*E1*E4^2 - 1.3201*E1*E4*(E5-1) - 0.092505*E1*(E5-1)^2 - 0.99636*E2^3 - 1.5818*E2^2*E4 - 2.6332*E2^2*(E5-1) - 1.5522*E2*E4^2 - 2.1611*E2*E4*(E5-1) - 2.3123*E2*(E5-1)^2 - 1.6751*E4^3 - 3.2083*E4^2*(E5-1) - 2.6873*E4*(E5-1)^2 - 3.2498*(E5-1)^3 - 0.071762*E1^2 + 0.054463*E1*E2 - 0.43302*E1*E4 + 0.053938*E1*(E5-1) + 0.93482*E2^2 + 0.27295*E2*E4 + 3.9166*E2*(E5-1) - 3.436*E4^2 - 7.4734*E4*(E5-1) + 0.77528*(E5-1)^2 + 2.233*E1 - 1.5702*E2 - 2.2243*E4 + 0.63052*(E5-1) + ubar(1);
    w(k) = -0.16412*E1^4 - 0.034221*E1^3*E2 - 0.21967*E1^3*E4 - 0.72915*E1^3*(E5-1) + 0.066252*E1^2*E2^2 + 0.4081*E1^2*E2*E4 + 0.62513*E1^2*E2*(E5-1) - 0.12956*E1^2*E4^2 - 0.1923*E1^2*E4*(E5-1) - 0.57963*E1^2*(E5-1)^2 - 0.075495*E1*E2^3 - 0.30878*E1*E2^2*E4 - 0.71119*E1*E2^2*(E5-1) - 0.12179*E1*E2*E4^2 - 0.028558*E1*E2*E4*(E5-1) + 0.13849*E1*E2*(E5-1)^2 - 0.20053*E1*E4^3 - 0.5691*E1*E4^2*(E5-1) - 0.094239*E1*E4*(E5-1)^2 - 0.74132*E1*(E5-1)^3 + 0.22033*E2^4 + 0.43021*E2^3*E4 + 0.51376*E2^3*(E5-1) + 0.27807*E2^2*E4^2 - 0.27766*E2^2*E4*(E5-1) - 0.18055*E2^2*(E5-1)^2 + 0.42377*E2*E4^3 + 0.45828*E2*E4^2*(E5-1) + 0.36612*E2*E4*(E5-1)^2 + 0.58741*E2*(E5-1)^3 + 0.065266*E4^4 - 0.10996*E4^3*(E5-1) - 0.18811*E4^2*(E5-1)^2 - 0.14792*E4*(E5-1)^3 - 0.36763*(E5-1)^4 - 0.14583*E1^3 - 0.068122*E1^2*E2 - 0.21675*E1^2*E4 - 0.75259*E1^2*(E5-1) - 0.07876*E1*E2^2 - 0.010202*E1*E2*E4 + 0.079325*E1*E2*(E5-1) - 0.20761*E1*E4^2 - 0.14919*E1*E4*(E5-1) - 0.09817*E1*(E5-1)^2 - 0.051686*E2^3 - 0.084017*E2^2*E4 - 0.5065*E2^2*(E5-1) - 0.10947*E2*E4^2 - 0.03038*E2*E4*(E5-1) - 0.037294*E2*(E5-1)^2 - 0.22573*E4^3 - 0.66462*E4^2*(E5-1) + 0.28994*E4*(E5-1)^2 - 0.45811*(E5-1)^3 - 0.1266*E1^2 + 0.043646*E1*E2 - 0.56*E1*E4 - 0.8092*E1*(E5-1) - 0.047322*E2^2 - 0.45767*E2*E4 - 0.38652*E2*(E5-1) - 1.1646*E4^2 - 3.5854*E4*(E5-1) - 1.7053*(E5-1)^2 - 0.19172*E1 + 1.2057*E2 + 1.9802*E4 + 4.1205*(E5-1) + ubar(2);
    %v(k) = K1*subs(Z)+ubar(1);
    %w(k) = K2*subs(Z)+ubar(2);

    X = [X x+dt*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
end

Error = Xr-X(:,1:end-1);
Error4 = sin(Error(3,:));
Error5 = cos(Error(3,:));

plot(X(1,1:end-1),X(2,1:end-1),'r')
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
plot(time,Error(1,:),time,Error(2,:),time,Error(3,:),time,Error4(:),time,Error5(:));
legend('e1','e2','e3','e4','e5')
