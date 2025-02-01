%% Parameters for Buck Converter Polynomial Voltage Control Example
% 
% This example shows how to control the output voltage of a buck 
% converter using a polynomial RST controller. The RST controller adjusts 
% the duty cycle. The input voltage is considered constant throughout the 
% simulation. A variable resistor provides the load for the system. The 
% total simulation time (t) is 0.25 seconds. At t = 0.15 seconds, the load 
% is changed. At t = 0.2 seconds, the voltage reference is changed from 
% 6V to 4V.

% Copyright 2017-2023 The MathWorks, Inc.

%% System Parameters
Vin = 12;     % Input voltage [V]
L   = 0.0037; % Inductance    [H]
C   = 1e-5;   % Capacitance   [F] 
R   = 0.01;   % Capacitor effective series resistance [Ohm]

%% Control Parameters
Ts  = 5e-6; % Fundamental sample time            [s]
Tsc = 1e-4; % Sample time for inner control loop [s]

Rpol = [1,-0.985076074360045,-0.0432890026766142,0.0283650770366593];
Spol = [0.0850919621380208,-0.121945676411019,0.0472050606201127];
Tpol = [0,0.00351848376928389,0.00683286257783044];
