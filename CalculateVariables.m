function [ variables, I, rates, fluxes, sensitivity, dIdp, voltage, dVdt, time, G ] = CalculateVariables( protocol, exp_ref, fitting_protocol )
%CALCULATEVARIABLES Summary of this function goes here
%   Detailed explanation goes here

% Identifies temperature for appropriate experiment
if nargin == 2
    fitting_protocol = 'sine_wave';
end

temperature = GetTemperature( exp_ref );

if strcmp(exp_ref,'16713110')==1
    temperature = 21.4;
end

% Identifies best fitting parameters to sine wave from Kylie's work
model = 'hh';
[chain,likelihood] = FindingBestFitsAfterMCMC( model, fitting_protocol, exp_ref);
[~,v]= max(likelihood);
params = chain(v,:);

%Import action potential protocol
cd ../Protocols
voltage=importdata([protocol '_protocol.mat']);
cd ..
cd Code
% Simulate model with parameters identified as providing bet fit to sine wave
[I, Y, sensitivity, dIdp, fluxes, time, G ] = SimulatingDataSens(protocol,params,voltage,temperature);
rates = CalculateRates( voltage, params );
variables = [ Y, 1 - sum( Y, 2 )  ];
dVdt = diff(voltage)/0.1;
dVdt = [ dVdt; dVdt(end) ];

end

