% Code to plot simulation results from BuckPolynomialControl
%% Plot Description:
%
% The plot below shows the requested and measured voltage for the
% test and the input voltage in the circuit.

% Copyright 2017-2023 The MathWorks, Inc.

% Generate simulation results if they don't exist
if ~exist('simlog_BuckPolynomialControl', 'var') || ...
        simlogNeedsUpdate(simlog_BuckPolynomialControl, 'BuckPolynomialControl') 
    sim('BuckPolynomialControl')
    % Model StopFcn callback adds a timestamp to the Simscape simulation data log
end


% Reuse figure if it exists, else create new figure
if ~exist('h1_BuckPolynomialControl', 'var') || ...
        ~isgraphics(h1_BuckPolynomialControl, 'figure')
    h1_BuckPolynomialControl = figure('Name', 'BuckPolynomialControl');
end
figure(h1_BuckPolynomialControl)
clf(h1_BuckPolynomialControl)

% Get simulation results
simlog_t = simlog_BuckPolynomialControl.Buck_Converter.p2.v.series.time;
simlog_Vout = simlog_BuckPolynomialControl.Sensing_Vout.Voltage_Sensor.V.series.values('V');
simlog_Vin = simlog_BuckPolynomialControl.Sensing_Vin.Voltage_Sensor.V.series.values('V');
simlog_vRef = logsout_BuckPolynomialControl.get('voltage_request');

% Plot results
simlog_handles(1) = subplot(2, 1, 1);
plot(simlog_t, simlog_Vout, 'LineWidth', 1)
hold on
plot(simlog_vRef.Values.Time, simlog_vRef.Values.Data, 'LineWidth', 1)
hold off
grid on
title('Output voltage')
ylabel('Voltage (V)')
legend({'Measured','Reference'},'Location','Best');

simlog_handles(2) = subplot(2, 1, 2);
plot(simlog_t, simlog_Vin, 'LineWidth', 1)
grid on
title('Input voltage')
ylabel('Voltage (V)')
xlabel('Time (s)')

linkaxes(simlog_handles, 'x')

% Remove temporary variables
clear simlog_t simlog_handles temp_colororder
clear simlog_vRef simlog_Vout simlog_Vin
