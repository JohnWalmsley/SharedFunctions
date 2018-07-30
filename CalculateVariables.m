function [ variables, I, rates, fluxes, sensitivity, dIdp, voltage, dVdt, time, G ] = CalculateVariables( protocol, exp_ref, fitting_protocol, model )
%CALCULATEVARIABLES Summary of this function goes here
%   Detailed explanation goes here

% Identifies temperature for appropriate experiment
if nargin < 3
    fitting_protocol = { 'sine_wave' };
end
if nargin < 4
    model = 'hh';
end

temperature = GetTemperature( exp_ref );

if ~isa( fitting_protocol, 'cell' )
    fitting_protocol={fitting_protocol};
end

% Identifies best fitting parameters to sine wave from Kylie's work
[chain,likelihood] = FindingBestFitsAfterMCMC( model, fitting_protocol, exp_ref);
[~,v]= max(likelihood);
params = chain(v,:);

%Import action potential protocol
cd ../Protocols
voltage=importdata([protocol '_protocol.mat']);
cd ..
cd Code

[ ~, model_type ] = modeldata( model );
% Simulate model with parameters identified as providing bet fit to sine wave
[I, Y, sensitivity, dIdp, fluxes, time, G ] = SimulatingDataSens(protocol,params,voltage,temperature, model_type);
rates = CalculateRates( voltage, params, model_type );
variables = [ Y, 1 - sum( Y, 2 )  ];
dVdt = diff(voltage)/0.1;
dVdt = [ dVdt; dVdt(end) ];

end

