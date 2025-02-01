

A = [1 0 2; 0 -1 1; 0 2 -1];
B = [1 1; 0 0; 0 2];
C = [1 0 0];
% x = [x1 x2 x3];

a = 5;
% [x1,x2,x3] = meshgrid(-a:1:a,-a:1:a,-a:1:a);
[x1,x3] = meshgrid(-a :1:a,-a:1:a);
x2=0;
dx1 = A(1,1).*x1 + A(1,2).*x2 + A(1,3).*x3;
dx2 = A(2,1).*x1 + A(2,2).*x2 + A(2,3).*x3;
dx3 = A(3,1).*x1 + A(3,2).*x2 + A(3,3).*x3;
quiver(x1,x3,dx1,dx3)

pvar x1 x2 x3;
x = [x1;x2;x3];
m=2;
n=3;
p=1;
prog = sosprogram(x);
[prog,P0]=sospolymatrixvar(prog,monomials(x,0),[n,n],'symmetric');
[prog,N0]=sospolymatrixvar(prog,monomials(x,0),[n,n],'symmetric');
[prog,Q0]=sospolymatrixvar(prog,monomials(x,0),[p,p],'symmetric');
[prog,S0]=sospolymatrixvar(prog,monomials(x,0),[p,m]); 
[prog,R0]=sospolymatrixvar(prog,monomials(x,0),[m,m],'symmetric');

cond=[A'*P0+P0*A+N0-C'*Q0*C P0*B-C'*S0; (P0*B-C'*S0)' -R0];
prog=sosineq(prog,-cond);
prog=sosineq(prog,P0);
prog=sosineq(prog,N0);
prog=sosineq(prog,R0);

sol=sossolve(prog);
P0=double(sosgetsol(sol,P0))
N0=double(sosgetsol(sol,N0))
Q0=double(sosgetsol(sol,Q0))
S0=double(sosgetsol(sol,S0))
R0=double(sosgetsol(sol,R0))

k = 0;
kmax = 20;
if (S0*inv(R0)*S0'-Q0) >= 0 
    k=kmax;
    K=-inv(R0)*S0'
    V0 = x'*P0*x
end
while k < kmax
prog = sosprogram(x);
[prog,P]=sospolymatrixvar(prog,monomials(x,0),[n,n],'symmetric');
[prog,N]=sospolymatrixvar(prog,monomials(x,0),[n,n],'symmetric');
[prog,Q]=sospolymatrixvar(prog,monomials(x,0),[p,p],'symmetric');
[prog,S]=sospolymatrixvar(prog,monomials(x,0),[p,m]); 
[prog,R]=sospolymatrixvar(prog,monomials(x,0),[m,m],'symmetric');

cond1=[A'*P+P*A+N-C'*Q*C P*B-C'*S; (P*B-C'*S)' -R];
cond2=R0-R;
cond3=S*inv(R0)*S0'+S0*inv(R0)*S'-2*S0*inv(R0)*S0'+Q0-Q;
prog=sosineq(prog,-cond1);
prog=sosineq(prog,cond2);
prog=sosineq(prog,cond3);
prog=sosineq(prog,P);
prog=sosineq(prog,N);
prog=sosineq(prog,R);

sol=sossolve(prog);
P=double(sosgetsol(sol,P))
N=double(sosgetsol(sol,N))
Q=double(sosgetsol(sol,Q))
S=double(sosgetsol(sol,S))
R=double(sosgetsol(sol,R))

Q0=Q;
S0=S;
R0=R;
k = k+1;
if S*inv(R)*S'-Q >= 0
    K=-inv(R)*S'
    break
end
end

% [x1,x2,x3] = meshgrid(-a:1:a,-a:1:a,-a:1:a);
[x1,x3] = meshgrid(-a:1:a,-a:1:a);
x2=0;
dx1 = A(1,1).*x1 + A(1,2).*x2 + A(1,3).*x3 + B(1,:)*K.*x1;
dx2 = A(2,1).*x1 + A(2,2).*x2 + A(2,3).*x3 + B(2,:)*K.*x1;
dx3 = A(3,1).*x1 + A(3,2).*x2 + A(3,3).*x3 + B(3,:)*K.*x1;
figure
quiver(x1,x3,dx1,dx3)

