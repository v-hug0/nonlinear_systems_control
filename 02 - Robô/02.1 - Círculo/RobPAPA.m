clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SUM-OF-SQUARES FEASIBILITY (OR OPTIMIZATION) PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initializations
pvar e1 e2 e4 e5
e_n=[e1;e2;e4;e5];
Ur = [0.2;0.2];
ubar = Ur;


A = [0 0 0 0;
     0 0 0 0;
     0 vr 0 wr;
     vr 0 wr 0];
B = [-1  e2;
      0 -e1;
      0 -e5-1;
      0  e4];


z = monomials([e1;e2;e4;e5+1],1);

n = length(e_n);
nz = length(nz);
m = size(z,2);

prog = sosprogram(e_n);
[prog,P] = sospolymatrixvar(prog,monomials(e_n,0),[nz,nz],'symmetric');
[prog,K] = sospolymatrixvar(prog,monomials(e_n,0:3),[m,nz]);
[prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
[prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)],'wscoeff');
% Defining the SOS constraints
g1 = e4^2+e5^2-1; % constraint from eq.(13)
gradV = [diff(V,e1);diff(V,e2);diff(V,e4);diff(V,e5)];
dV = 
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
[prog,l1] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
[prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
[prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)], 'wscoeff');
% Defining the SOS constraints
g1 = e4^2+e5^2-1; % constraint from eq.(13)
gradV = [diff(V,e1);diff(V,e2);diff(V,e4);diff(V,e5)];
dV = -gradV'*(f+g*(u+ubar))-s1*g1-l2;
prog = sosineq(prog,dV);
prog = sosineq(prog,V-l1);
% Second solution: (V,l1,l2,s1)
sol = sossolve(prog);
V = sosgetsol(sol,V);
l1 = sosgetsol(sol,l1);
l2 = sosgetsol(sol,l2);
s1 = sosgetsol(sol,s1);


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
wr = 0.3;
vr = wr;
r = 1;

thetar = pi/2+wr*t;
xr = r*cos(wr*t);
yr = r*sin(wr*t);
Xr=[xr;yr;thetar];

figure
plot(Xr(1,1),Xr(2,1),'o','MarkerFaceColor',	"#4DBEEE");
hold on
plot(Xr(1,:),Xr(2,:),'LineWidth',2,'Color',	"#4DBEEE");
hold on
%%

h = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          ROBOT REFERENCE TRACKING SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [0;0;pi/4];
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

    v(k) = K1*subs(Z)+ubar(1);
    w(k) = K2*subs(Z)+ubar(2);

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
