function jacobian = CalculateJacobianMatrix( v, params )
%CALCULATEJACOBIANMATRIX Summary of this function goes here
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

r1 = [ -( k12 + k14 ) - k41,  k21 - k41,  -k41 ];
r2 = [ k12, -( k23 + k21 ), k32 ];
r3 = [ -k43, k23-k43, -( k34 + k32 ) - k43 ];

jacobian = [r1; r2; r3 ];

end

