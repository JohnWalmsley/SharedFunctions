function [ exp_data_corrected, Ileak ] = CorrectLeakCurrent( V, exp_data, leak_idxs )
%CORRECTLEAKCURRENT Corrects data for leak
%   This code calculates the leak current based on the leak step at the
%   beginning of each protocol. 

% We require that the current is zero at the beginning of the trace to
% match the data!

v_rest = mean( V( 1 : leak_idxs( 1 ) - 100 ) );
v_step = mean( V( leak_idxs) );

Ileak_rest = mean( exp_data( 1 : leak_idxs( 1 ) - 100 ) );
Ileak_step = mean( exp_data( leak_idxs ) );

% I leak = Gleak * V - const 

params = [ v_rest -1; v_step -1 ] \ [ Ileak_rest ; Ileak_step ];

Ileak = params( 1 ) * V - params( 2 );
exp_data_corrected = exp_data - Ileak;

end

