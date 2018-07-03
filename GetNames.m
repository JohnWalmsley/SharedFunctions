function Names = GetNames( conf )
%GETNAMES Summary of this function goes here
%   Detailed explanation goes here

variable_names = { 'y1', 'y2', 'y3', 'y4' };
I_names = { 'I' };
rates_names = { 'k12', 'k14', 'k21', 'k23', 'k32', 'k34', 'k41', 'k43' };
fluxes_names = { 'k12*y1', 'k14*y1', 'k21*y2', 'k23*y2', 'k32*y3', 'k34*y3', 'k41*y4', 'k43*y4' };
sensitivity_names = {};
for var = 1 : 3
    for par = 0 : 7
        sensitivity_names{ 8*(var-1) + par+1 } = [ 'dy' num2str( var ) '/dp' num2str( par ) ];
    end
end
dIdp_names = { 'dIdp0','dIdp1','dIdp2','dIdp3','dIdp4','dIdp5','dIdp6','dIdp7' };
voltage_names = {'V'};
dVdt_names = {'dV/dt'};
const_names = {'const'};
time_names = {'time'};

Names = {};

if any(strcmp( conf, 'const' ) )
    Names = [ Names const_names ];
end
if any(strcmp( conf, 'variables' ) )
    Names = [ Names variable_names ];
end
if any(strcmp( conf, 'I' ) )
   Names = [ Names I_names ];
end
if any(strcmp( conf, 'rates' ) )
   Names = [ Names rates_names ];
end
if any(strcmp( conf, 'fluxes' ) )
   Names = [ Names fluxes_names ];
end
if any(strcmp( conf, 'sensitivity' ) )
   Names = [ Names sensitivity_names ];
end
if any(strcmp( conf, 'dIdp' ) )
   Names = [ Names dIdp_names ];
end
if any(strcmp( conf, 'voltage' ) )
   Names = [ Names voltage_names ];
end
if any(strcmp( conf, 'dVdt' ) )
   Names = [ Names dVdt_names ];
end
if any(strcmp( conf, 'time' ) )
   Names = [ Names time_names ];
end

end

