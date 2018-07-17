function rhsVectors = CalculateSensitivityRhsVectors( v, y, params, model_type )
%UNTITLED Summary of this function goes here
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

if model_type == 35
    
    k12 = k43;
    k21 = k34;
    
    k41 = k32;
    k14 = k23;
    
    y1 = y(1);  y2 = y(2);  y3 = y(3);
    y4 = 1 - y1 - y2 - y3;
    
    c1 = [ -exp(P1*v)*y1; exp(P1*v)*y1; exp(P1*v)*y4 ];
    c2 = [ -v*k12*y1; v*k12*y1; v*k43*y4 ];
    c3 = [ exp(-P3*v)*y2; -exp(-P3*v)*y2; -exp(-P3*v)*y3 ];
    c4 = [ -v*k21*y2; v*k21*y2;  v*k34*y3 ];
    c5 = [  exp(P5*v)*y4; exp(P5*v)*y3; -exp(P5*v)*y3; ];
    c6 = [ v*k41*y4; v*k32*y3; -v*k32*y3 ];
    c7 = [ exp(-P7*v)*y1; -exp(-P7*v)*y2; exp(-P7*v)*y2 ];
    c8 = [ v*k14*y1; v*k23*y2; -v*k23*y2 ];
    
    rhsVectors = [c1 c2 c3 c4 c5 c6 c7 c8 ];

elseif any( model_type == [ 36 37 38 39 40 ] )
    
    m = y( 1 );
    h = y( 2 );
    
    k1 = k43;
    k2 = k34;
    k3 = k32;
    k4 = k23;
    
    c1 = [ ( 1 - m ) * exp( P1 * v ); 0 ];
    c2 = [ ( 1 - m ) * v * k1; 0 ];
    c3 = [ -m * exp( -P3 * v ); 0 ];
    c4 = [ m * v * k2; 0 ];
    c5 = [ 0; ( 1 - h ) * exp( P5 * v ) ]; 
    c6 = [ 0; ( 1 - h ) * v * k3 ];
    c7 = [ 0; -h * exp( -P7 * v ) ];
    c8 = [ 0;  h * v * k4 ];

    rhsVectors = [c1 c2 c3 c4 c5 c6 c7 c8 ];
end

