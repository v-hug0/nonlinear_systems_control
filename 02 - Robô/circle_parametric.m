clear
clc

a = ;
syms t
% Compute x and y using parametric equations
x = a*cos(t);
y = a*sin(t);

% Compute derivatives of x and y with respect to time
dx = diff(x,t);
dy = diff(y,t);
v = sqrt(dx^2 + dy^2);

% Compute angular velocity (rate of change of phi with respect to time)
d2x = diff(dx,t);
d2y = diff(dy,t);
%theta = atan2(dy,dx);
K = 1/R;
dtheta = K*v;
w = (dx*d2y-dy*d2x)/(dx^2+dy^2);


%%

num_points = 1000;  % number of points
k = 1;

Ts = 2*pi;
dt=0.01;

a = 1;
% Compute theta values
t = linspace(0, 2*pi*k, num_points);  % linear theta values from 0 to 2*pi;
xr = a*cos(t)+1;
yr = a*sin(t)+1;
theta_r = atan2(-cos(t),sin(t))+k*pi;
Xr = [xr;yr;theta_r];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Simulacao do robo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xo = [0;0;0];
%xo = [x(1)-1;y(1)-1;theta(1)];
x = xo;
X = x;
v = [];
w = [];
ubar = [1;1];

for k = 1:(num_points-1)    % 1 a tempo de simula��o/ tempo de integra��o  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(1);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);

    %v(k) = 0.022679*E1^4 + 0.023428*E1^3*E2 + 0.060221*E1^3*E4 - 0.014796*E1^3*(E5-1) - 0.07072*E1^2*E2^2 - 0.048277*E1^2*E2*E4 - 0.0029425*E1^2*E2*(E5-1) + 0.026879*E1^2*E4^2 - 0.016501*E1^2*E4*(E5-1) + 0.0049379*E1^2*(E5-1)^2 - 0.0090439*E1*E2^3 - 0.013502*E1*E2^2*E4 + 0.031065*E1*E2^2*(E5-1) - 0.012775*E1*E2*E4^2 + 0.0098654*E1*E2*E4*(E5-1) - 0.0091849*E1*E2*(E5-1)^2 + 0.0033585*E1*E4^3 - 0.004009*E1*E4^2*(E5-1) - 0.00030846*E1*E4*(E5-1)^2 - 2.8894E-05*E1*(E5-1)^3 - 6.9073E-05*E2^4 - 0.00032991*E2^3*E4 + 0.0025334*E2^3*(E5-1) - 0.001316*E2^2*E4^2 + 0.0020281*E2^2*E4*(E5-1) - 0.0034677*E2^2*(E5-1)^2 - 0.00045076*E2*E4^3 + 0.001667*E2*E4^2*(E5-1) - 0.0010074*E2*E4*(E5-1)^2 + 0.0018205*E2*(E5-1)^3 + 8.2562E-05*E4^4 - 0.00014535*E4^3*(E5-1) - 0.00023753*E4^2*(E5-1)^2 + 0.00016903*E4*(E5-1)^3 - 0.00021834*(E5-1)^4 + 0.63377*E1^3 + 0.11797*E1^2*E2 - 0.17067*E1^2*E4 - 0.032735*E1^2*(E5-1) + 0.59389*E1*E2^2 - 0.12065*E1*E2*E4 + 0.094889*E1*E2*(E5-1) + 0.2795*E1*E4^2 + 0.060938*E1*E4*(E5-1) + 0.36324*E1*(E5-1)^2 - 0.067703*E2^3 + 0.17016*E2^2*E4 + 0.033946*E2^2*(E5-1) + 0.025644*E2*E4^2 + 0.049392*E2*E4*(E5-1) + 0.019576*E2*(E5-1)^2 - 0.083759*E4^3 - 0.094921*E4^2*(E5-1) - 0.069397*E4*(E5-1)^2 - 0.0077285*(E5-1)^3 + 0.94075*E1^2 + 0.021785*E1*E2 + 0.046243*E1*E4 + 0.049291*E1*(E5-1) + 0.28182*E2^2 + 0.22786*E2*E4 - 0.0076298*E2*(E5-1) + 0.2202*E4^2 - 0.047299*E4*(E5-1) - 0.097279*(E5-1)^2 + 0.69754*E1 - 0.053633*E2 - 0.12827*E4 + 0.18276*(E5-1) + ubar(1);
    %w(k) = -0.08465*E1^4 + 0.10397*E1^3*E2 - 0.043461*E1^3*E4 + 0.031665*E1^3*(E5-1) + 0.013109*E1^2*E2^2 + 0.022022*E1^2*E2*E4 - 0.040742*E1^2*E2*(E5-1) - 0.01057*E1^2*E4^2 + 0.011835*E1^2*E4*(E5-1) - 0.0053815*E1^2*(E5-1)^2 + 0.0078426*E1*E2^3 - 0.00034521*E1*E2^2*E4 - 0.0038745*E1*E2^2*(E5-1) + 0.0073734*E1*E2*E4^2 - 0.0048454*E1*E2*E4*(E5-1) + 0.013511*E1*E2*(E5-1)^2 - 0.0018602*E1*E4^3 + 0.0014257*E1*E4^2*(E5-1) - 0.0024962*E1*E4*(E5-1)^2 - 0.00013949*E1*(E5-1)^3 + 0.00087367*E2^4 + 0.00026746*E2^3*E4 - 0.0011089*E2^3*(E5-1) + 0.00069883*E2^2*E4^2 + 1.7594E-05*E2^2*E4*(E5-1) + 0.0015514*E2^2*(E5-1)^2 + 0.00033788*E2*E4^3 - 0.00099626*E2*E4^2*(E5-1) + 0.00062091*E2*E4*(E5-1)^2 - 0.0014278*E2*(E5-1)^3 - 9.9931E-05*E4^4 + 7.9344E-05*E4^3*(E5-1) + 0.00013545*E4^2*(E5-1)^2 + 8.4911E-05*E4*(E5-1)^3 + 0.00035443*(E5-1)^4 + 0.14598*E1^3 - 0.66552*E1^2*E2 + 0.48661*E1^2*E4 - 0.19094*E1^2*(E5-1) - 0.2131*E1*E2^2 - 0.10649*E1*E2*E4 - 0.17653*E1*E2*(E5-1) + 0.15813*E1*E4^2 + 0.11582*E1*E4*(E5-1) + 0.23755*E1*(E5-1)^2 + 0.92623*E2^3 - 0.62134*E2^2*E4 + 0.33491*E2^2*(E5-1) - 0.31875*E2*E4^2 - 0.17727*E2*E4*(E5-1) - 0.43811*E2*(E5-1)^2 + 0.32168*E4^3 + 0.017083*E4^2*(E5-1) - 0.081435*E4*(E5-1)^2 - 0.039005*(E5-1)^3 - 0.10656*E1^2 + 0.23632*E1*E2 - 0.59375*E1*E4 + 0.16229*E1*(E5-1) + 0.014501*E2^2 - 0.0059766*E2*E4 + 0.12526*E2*(E5-1) + 0.032654*E4^2 - 0.067666*E4*(E5-1) + 0.11975*(E5-1)^2 - 0.18687*E1 + 1.3487*E2 + 0.9631*E4 + 0.08193*(E5-1) + ubar(2);

    v(k) = -0.0050973*E1^4 + 0.038267*E1^3*E2 + 0.0071696*E1^3*E4 + 0.011913*E1^3*(E5-1) - 0.013063*E1^2*E2^2 - 0.059192*E1^2*E2*E4 - 0.011194*E1^2*E2*(E5-1) - 0.010027*E1^2*E4^2 - 0.012052*E1^2*E4*(E5-1) - 0.011189*E1^2*(E5-1)^2 - 0.025969*E1*E2^3 + 0.022535*E1*E2^2*E4 - 0.026544*E1*E2^2*(E5-1) + 0.032729*E1*E2*E4^2 + 0.0043205*E1*E2*E4*(E5-1) - 0.0073466*E1*E2*(E5-1)^2 + 0.0075101*E1*E4^3 + 0.0051714*E1*E4^2*(E5-1) + 0.010027*E1*E4*(E5-1)^2 - 0.000175*E1*(E5-1)^3 - 2.4068E-05*E2^4 + 0.0090417*E2^3*E4 - 0.0013576*E2^3*(E5-1) - 0.00059123*E2^2*E4^2 + 0.0049751*E2^2*E4*(E5-1) + 0.0056068*E2^2*(E5-1)^2 - 0.006802*E2*E4^3 + 0.0062368*E2*E4^2*(E5-1) - 0.0029803*E2*E4*(E5-1)^2 + 0.005797*E2*(E5-1)^3 - 0.0021596*E4^4 + 0.0016337*E4^3*(E5-1) - 0.0027679*E4^2*(E5-1)^2 + 0.0013411*E4*(E5-1)^3 + 0.0020022*(E5-1)^4 + 2.1959*E1^3 + 0.22746*E1^2*E2 + 1.2948*E1^2*E4 - 0.06067*E1^2*(E5-1) + 0.38652*E1*E2^2 + 0.0045471*E1*E2*E4 - 0.41956*E1*E2*(E5-1) + 1.8203*E1*E4^2 + 0.30167*E1*E4*(E5-1) + 1.586*E1*(E5-1)^2 + 0.019691*E2^3 - 0.14687*E2^2*E4 - 0.23582*E2^2*(E5-1) + 0.30674*E2*E4^2 - 0.13482*E2*E4*(E5-1) + 0.2797*E2*(E5-1)^2 + 0.049388*E4^3 + 0.37237*E4^2*(E5-1) + 0.52972*E4*(E5-1)^2 + 0.11415*(E5-1)^3 + 1.6958*E1^2 - 0.097993*E1*E2 + 0.42105*E1*E4 - 0.3504*E1*(E5-1) - 0.62425*E2^2 - 0.19144*E2*E4 + 0.13537*E2*(E5-1) - 0.078572*E4^2 - 0.11852*E4*(E5-1) - 0.10201*(E5-1)^2 + 1.6464*E1 + 0.15451*E2 + 0.12518*E4 + 0.56361*(E5-1) + ubar(1);
    w(k) = -0.085834*E1^4 + 0.042031*E1^3*E2 + 0.18548*E1^3*E4 + 0.026801*E1^3*(E5-1) + 0.054055*E1^2*E2^2 - 0.084631*E1^2*E2*E4 + 0.036047*E1^2*E2*(E5-1) - 0.16932*E1^2*E4^2 - 0.0083642*E1^2*E4*(E5-1) - 0.0033644*E1^2*(E5-1)^2 + 0.0020219*E1*E2^3 - 0.047399*E1*E2^2*E4 + 0.0034378*E1*E2^2*(E5-1) + 0.059386*E1*E2*E4^2 - 0.033286*E1*E2*E4*(E5-1) + 0.008371*E1*E2*(E5-1)^2 + 0.074279*E1*E4^3 - 0.010661*E1*E4^2*(E5-1) + 0.0023638*E1*E4*(E5-1)^2 + 0.0019696*E1*(E5-1)^3 + 0.0051378*E2^4 - 0.0015368*E2^3*E4 + 0.0081664*E2^3*(E5-1) + 0.02306*E2^2*E4^2 - 0.011125*E2^2*E4*(E5-1) + 0.016481*E2^2*(E5-1)^2 - 0.015409*E2*E4^3 + 0.020115*E2*E4^2*(E5-1) - 0.01399*E2*E4*(E5-1)^2 + 0.011823*E2*(E5-1)^3 - 0.013017*E4^4 + 0.0098133*E4^3*(E5-1) - 0.00095033*E4^2*(E5-1)^2 + 0.0011707*E4*(E5-1)^3 + 0.003789*(E5-1)^4 - 0.45406*E1^3 + 0.25898*E1^2*E2 - 0.1644*E1^2*E4 - 0.45515*E1^2*(E5-1) + 0.13519*E1*E2^2 + 0.70012*E1*E2*E4 - 0.21406*E1*E2*(E5-1) - 0.72342*E1*E4^2 + 0.014562*E1*E4*(E5-1) - 0.61156*E1*(E5-1)^2 + 4.1129*E2^3 + 0.30872*E2^2*E4 - 1.6564*E2^2*(E5-1) + 0.046438*E2*E4^2 + 0.076675*E2*E4*(E5-1) - 0.02696*E2*(E5-1)^2 + 0.72704*E4^3 + 0.33709*E4^2*(E5-1) + 0.52808*E4*(E5-1)^2 + 0.19449*(E5-1)^3 - 0.54554*E1^2 + 0.0077779*E1*E2 - 0.094133*E1*E4 + 0.15673*E1*(E5-1) - 0.55987*E2^2 - 0.12359*E2*E4 + 0.30433*E2*(E5-1) + 0.13966*E4^2 + 0.28297*E4*(E5-1) + 0.14675*(E5-1)^2 + 0.057877*E1 + 2.4576*E2 + 1.4558*E4 - 0.10733*(E5-1) + ubar(2);

    X = [X x+dt*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
    
    Erro(:,k) = Xr(:,k+1)- X(:,k+1);
   
end

plot(xr,yr);
hold on
plot(X(1,:),X(2,:),'r')
legend('rungeKutta','trajetoria do robo')
xlabel('x[m]')
ylabel('y[m]')

figure
v = (cos(t).^2 + sin(t).^2).^(1/2);
%v = ((4*cos(t).^2)/25 + (4*sin(t).^2)/25).^(1/2)
subplot(2,1,2)
plot(t,v)
w = zeros(1,length(t))+w;
subplot(2,1,2)
plot(t,w)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem de v e w 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
ne = length(Erro);
tempo = Ts*(0:ne-1);
subplot(2,1,1)
plot(tempo,v,'b')
legend('v')

hold on
subplot(2,1,2)
plot (tempo,w,'r')
legend('w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem do erro x y e theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(tempo,Erro(1,:),tempo,Erro(2,:),tempo,Erro(3,:));
legend('x','y','theta')

 EE1 = Erro(1,:);
 EE2 = Erro(2,:);
 EE3 = Erro(3,:);