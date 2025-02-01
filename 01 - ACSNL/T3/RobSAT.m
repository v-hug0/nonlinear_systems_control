clear;close all;clc

% Inicializações
pvar e1 e2 e4 e5
e = [e1;e2;e4;e5];
ebar = [0;0;0;1];
Ur = [0.7;0];
ubar = Ur;
% Input-affine Robot tracking error dynamics
f = [Ur(1)*e5+e2*ubar(2); Ur(1)*e4-e1*ubar(2); Ur(2)*e5; 0];
g = [-1 0.2*e2; 0 -0.2*e1; 0 -0.2*e5-0.2*1; 0 0.2*e4];
% Vector of monomials chosen
Z = e;
% Linear-like Robot tracking error dynamics 
A =  [0 ubar(2) 0 Ur(1);
      -ubar(2) 0 Ur(1) 0;
      0 0 0 Ur(2);
      0 0 0 0];
B = g;
% Dimensions
n = length(e);
m = size(g,2);
nz = length(Z);
% definicao de X = E(Sx)
Sx = eye(nz)*3;
epsi = 1e-2;
% -- SOS PROGRAM
prog = sosprogram(e);
[prog,Q] = sospolymatrixvar(prog,monomials(e,0),[nz,nz],'symmetric');
[prog,K] = sospolymatrixvar(prog,monomials(e,0:1),[m,nz]);
[prog,T] = sospolymatrixvar(prog,monomials(e,0:1),[m,nz]);
gx = 1-Z'*Sx*Z; %elipsoide of X
geq = e4^2+e5^2-1; % restrição de igualdade
for i=1:nz
    for j=1:n
    M(i,j) = diff(Z(i),e(j)); % matriz M (derivadas parciais de Z(x))
    end
end
% -- Theorem 1 - Eq.(09)
% vertΩ set of θ
th1 = [0;0];
th2 = [0;1];
th3 = [1;0];
th4 = [1;1];
vert_omega = [th1 th2 th3 th4];
n_vert = 2^m; % n° vertices of omega
for i=1:n_vert
    [prog,Sr] = sospolymatrixvar(prog,monomials(e,1:2),[nz,nz],'symmetric');
    [prog,Seq] = sospolymatrixvar(prog,monomials(e,0:1),[nz,nz],'symmetric');
    TH = eye(m)*vert_omega(i);
    F1 = -(M*A*Q+M*B*(TH*K+(eye(m)-TH)*T) + Q*A'*M' + (K'*TH'+T'*(eye(m)-TH'))*B'*M') - Sr*gx - Seq*geq;
    S1 = F1 - epsi*eye(nz);
    prog = sosineq(prog,S1);
end
% Theorem 1 - Eq.(10)
for j=1:m
    [prog,Sh] = sospolymatrixvar(prog,monomials(e,1:2),[nz,nz],'symmetric');
    [prog,Seq] = sospolymatrixvar(prog,monomials(e,0:1),[nz,nz],'symmetric');
    F2 = Q-Sh*gx-Seq*geq;
    S2 = [1 T(j); T(j)' F2(j)];
    prog = sosineq(prog,S2);
end
% -- Teorema 1 - Eq.(11)
[prog,Seq] = sospolymatrixvar(prog,monomials(e,0:1),[nz,nz],'symmetric');
S3 = inv(Sx)-Q-Seq*geq;
prog = sosineq(prog,S3);
% -- Maximização do elipsoide de V
obj = trace(Q);
prog = sossetobj(prog,-obj);
% Solution
sol = sossolve(prog);
Q = double(sosgetsol(sol,Q));
K = sosgetsol(sol,K);
T = sosgetsol(sol,T);
% Lyapunov and Controllers
V = Z'*inv(Q)*Z;
F = K*inv(Q);
H = T*inv(Q);
u = F*Z;
v = H*Z;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           REFERENCE TRAJECTORY GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ur = [0.7;0];
ubar = Ur;

dt = 0.1;          % sampling time

Tsim = 60;   %90      % tempo de simulacao
h = 0.1;           % integrating time


xr = [2;3;0];
Xr = xr;         %Euler


for k = 1:round(Tsim/h)  % 1 a tempo de simulacao/ tempo de integracao   
    vr(k) = 0.7;
    wr(k) = 0;
    Xr = [Xr xr+h*([vr(k)*cos(xr(3));vr(k)*sin(xr(3));wr(k)])];  %Euler
    xr = Xr(:,k+1);
end
Ur = [vr;wr];
hold on
plot(Xr(1,:),Xr(2,:),'')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          ROBOT REFERENCE TRACKING SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [0;-5;1.1];
x = xo;
X = x;
v = [];
w = [];

for k = 1:(round(Tsim/h))    % 1 a tempo de simulacao/tempo de integracao  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(1);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);

    v(k) = 8.9174E-05*E1^2 + 6.6914E-05*E1*E2 - 0.00024889*E1*E4 + 0.0008632*E1*(E5-1) + 0.14416*E2^2 + 0.54212*E2*E4 - 9.5307E-06*E2*(E5-1) + 0.011349*E4^2 + 5.7297E-06*E4*(E5-1) + 0.00033398*(E5-1)^2 + 1.7733*E1 + 0.0012801*E2 + 0.00012541*E4 + 0.65664*(E5-1) + ubar(1);
    v(k) = min(1,max(-1,v(k)));
    w(k) = 0.00032696*E1^2 - 0.52926*E1*E2 - 0.1817*E1*E4 + 0.00011986*E1*(E5-1) - 0.00036744*E2^2 - 0.00024061*E2*E4 + 0.29971*E2*(E5-1) + 7.3019E-05*E4^2 - 1.0948*E4*(E5-1) - 3.299E-05*(E5-1)^2 + 4.2933E-05*E1 + 0.82341*E2 + 2.7667*E4 - 8.1838E-05*(E5-1) + ubar(2);
    w(k) = min(0.2,max(-0.2,w(k)));

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




