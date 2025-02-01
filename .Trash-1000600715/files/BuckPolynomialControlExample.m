%% Buck Converter Polynomial Voltage Control
% 
% This example shows how to control the output voltage of a buck 
% converter using a polynomial RST controller. The RST controller adjusts 
% the duty cycle. The input voltage is considered constant throughout the 
% simulation. A variable resistor provides the load for the system. The 
% total simulation time (t) is 0.25 seconds. At t = 0.15 seconds, the load 
% is changed. At t = 0.2 seconds, the voltage reference is changed from 
% 6V to 4V.
% 

% Copyright 2017-2023 The MathWorks, Inc.


%% Model

open_system('BuckPolynomialControl')

set_param(find_system('BuckPolynomialControl','FindAll', 'on','type','annotation','Tag','ModelFeatures'),'Interpreter','off')

%% Simulation Results from Simscape Logging
%%
%
% The plot below shows the requested and measured voltage for the
% test and the input voltage in the circuit.
%


BuckPolynomialControlPlotVoltage;

%%

