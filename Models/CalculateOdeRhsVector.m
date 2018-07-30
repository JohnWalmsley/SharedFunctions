function ode_vector = CalculateOdeRhsVector( v, params )
%CALCULATEODERHSVECTOR Summary of this function goes here
%   Detailed explanation goes here

P0 = params(2);
P1 = params(3);
P4 = params(6);
P5 = params(7);

k32 = P4*exp(P5*v);
k43 = P0*exp(P1*v);    
k41 = k32;

ode_vector = [ k41; 0; k43 ];

end

