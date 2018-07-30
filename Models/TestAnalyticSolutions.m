function TestAnalyticSolutions

% This function tests the analytic solutions to the model ODEs and the
% sensitvity equations against their analytic solutions.
addpath ..
load TestParams % Best result from Kylie's MCMC

% Set voltages for testing
voltage_vector = -120 : 20 : 40;

fig = figure;
I = eye(3);
for v = 1 : length(voltage_vector)
    
    % Test solution of model equations
    IC = [0.0,0.0,0.0];
    vv = v * ones( 1, 10000 );
    ProtocolLength = length( vv )/10;
    params = [ 1000 * voltage_vector(v), P ];
    time = (0:0.1:ProtocolLength)';
    [y_sim,~,flux_sim] = MexHHSens( time' , IC, ...
                                    params, zeros( 1,(length(params)-2)*length(IC)));  
    % Calculate analytic solution
    % Initial condition is [ 0 0 0 ] - does not feature in solution
    B_V = CalculateOdeMatrix( voltage_vector(v), params );
    l_V = CalculateOdeRhsVector( voltage_vector(v), params );
    % Calculate the analytic solution
    y_star = -inv(B_V)*l_V;
    y_analytic = zeros( length( time ), 3 );
    for i = 1:length(time)
        t =time( i );
        y_temp = ( I - expm( t * B_V ) ) * y_star;        
        y_analytic( i, : ) = y_temp'; % first sensitivities for y1, then y2, then y3 as in code
    end
    
    subplot( 1, 3, 1 )
    hold on
    plot( time(1:end-1), y_sim, 'r')
    plot( time, y_analytic, 'k--')
    hold off
    
    rates = CalculateRates( voltage_vector(v), params );

    y4_analytic = 1 - sum( y_analytic, 2 );
    
    % From MexHHSens.c:
    %   flux[0] = k12*y1; //IC->I
    %   flux[1] = k14*y1; //IC->C
    %   flux[2] = k21*y2; // I->IC 
    %   flux[3] = k23*y2; // I->O
    %   flux[4] = k32*y3; // O->I
    %   flux[5] = k34*y3; // O->C
    %   flux[6] = k41*y4; // C->IC
    %   flux[7] = k43*y4; // C->O
    flux_analytic = zeros( length( time ) , 8 );
    Y = [y_analytic, y4_analytic ];
    flux_analytic( :, [ 1 3 5 7 ] ) = bsxfun( @times, rates( [ 1 3 5 7 ] ) , Y );
    flux_analytic( :, [ 2 4 6 8 ] ) = bsxfun( @times, rates( [ 2 4 6 8 ] ) , Y );
      
    subplot( 1, 3, 2 )
    hold on
    plot( time(1:end-1), flux_sim, 'r') % Note - does not return the initial flux
    plot( time, flux_analytic, 'k--')
    hold off
    
    % Test solution of Sensitivity equations
    % Model initial conditions are the analytic steady state
    % Sensitivity initial condition is zero
    
    [ ~, s_sim , ~ ] = MexHHSens( time', y_star', ...
                            params, zeros( 1,(length(params)-2)*length(y_star)));
    % add in the initial conditions
    
    A_V = CalculateJacobianMatrix( voltage_vector(v), params );
    k_V = CalculateSensitivityRhsVectors( voltage_vector(v), y_star', params, 35 );
    s_star = -inv(A_V)*k_V;
    s_analytic = zeros( length( time ), 24 );
    for i = 1:length( time )
        t = time( i );
        s_analytic_mx = ( I - expm( t * A_V ) ) * s_star; 
        s_temp = s_analytic_mx';
        s_analytic( i, : ) = s_temp(:)'; % first sensitivities for y1, then y2, then y3 as in code
    end
    subplot( 1, 3, 3 )
    hold on
    plot( time(1:end-1), s_sim, 'r')
    plot( time, s_analytic, 'k--')
    hold off
    
end

subplot( 1,3,1)
hold on
xlabel( 'Time (ms)' )
ylabel( 'y (-)' ) 
hold off

subplot( 1,3,2)
hold on
xlabel( 'Time (ms)' )
ylabel( 'flux (-/ms)' )
hold off

subplot( 1,3,3)
hold on
xlabel( 'Time (ms)' )
ylabel( 'Sens (depends)' )
hold off