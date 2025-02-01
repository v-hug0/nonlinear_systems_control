%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            COMPLEIB - AC11                              %
%                       REALIMENTAÇÃO DE ESTADOS                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear all

%------------------------- MODELO DA PLANTA ------------------------------%
A = [-1.341 0.9933 0 -0.1689 -0.2518;
    43.223 -0.8693 0 -17.251 -1.5766;
    1.341 0.0067 0 0.1689 0.2518;
    0 0 0 -20 0;
    0 0 0 0 -20];
B = [0 0; 0 0; 0 0; 20 0; 0 20];

disp('Autovalores de A:')
autoMA = eig(A)
disp('Instável em malha aberta!')


%----------------------- DECLARAÇÃO DAS VARIÁVEIS ------------------------%
W = sdpvar(5,5,'symmetric','real');
L = sdpvar(2,5,'full','real');


%----------------------- CONDIÇÕES DE ESTABILIDADE -----------------------%
CE = [];   
CE = CE + ( W - (10^-5)*eye(5) >= 0 );
CE = CE + ( (W'*A'+A*W+L'*B'+B*L) - (10^-5)*eye(5) <= 0 );


%--------------- RESOLVENDO A LMI COM TODAS AS RESTRIÇÕES ----------------%
options = sdpsettings('solver','sedumi','verbose',0);
optimize(CE);


%---------------------- EXIBINDO A SOLUÇÃO OBTIDA ------------------------%
W = value(W);
L = value(L);

P = inv(W)
K = L*P

disp('Autovalores de A+BK:')
autoMF = eig(A+B*K)
disp('Estável em malha fechada!')


