function rates = CalculateRates( v, params )
%CALCULATERATES Summary of this function goes here
%   Detailed explanation goes here

P0 = params(2);
P1 = params(3);
P2 = params(4);
P3 = params(5);
P4 = params(6);
P5 = params(7);
P6 = params(8);
P7 = params(9);
  
k32 = P4*exp(P5*v);
k23 = P6*exp(-P7*v);
    
k43 = P0*exp(P1*v);
k34 = P2*exp(-P3*v);
    
k12 = k43;
k21 = k34;
    
k41 = k32;
k14 = k23;

rates = [ k12 k14 k21 k23 k32 k34 k41 k43 ];

end

