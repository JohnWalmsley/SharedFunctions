function CalculateImposedVoltage( protocol, tau )

if nargin == 1
    tau = 0.26;
end

load ( [ '../Protocols/' protocol '_protocol.mat' ] );

time = 0 : 0.1: length( T )/10 - 0.1;
vprotocol = spline( time, T );

figure; plot( time, T, 'k' );
options = odeset( 'MaxStep', 10 );
[ ~, T ] = ode23t(@(t,V) oderhs( t, V, vprotocol, tau ), time, T(1), options); % Should output at same time points as original file

save( [ '../Protocols/' protocol '_protocol_lag.mat' ], 'T' ); % Note - Kylie's files save voltage as T.

hold on
plot( time, T, 'r' );
plot( time, ppval( time, vprotocol ), 'b' );
hold off
xlabel( 'Time (ms) ' )
ylabel( 'Voltage (ms)' )

end

function dVdt = oderhs( t, V, vprotocol, tau )

dVdt = ( ppval( t, vprotocol ) - V ) / tau;

end