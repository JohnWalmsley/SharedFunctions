function idxs = DetermineNoSpikeIdxs( data, thresh )

if nargin == 1
    thresh = 1;
end

diff_data = diff( data );

data_over_idx = find( abs( diff_data ) > thresh ); 

delete_idx = [];
for i = 2 : length( data_over_idx )
    if data_over_idx( i ) < data_over_idx( i-1) + 50
        delete_idx = [ delete_idx i];
    end
end

data_over_idx( delete_idx ) = [];
data_over_idx = data_over_idx-1;

idxs = zeros( length( data_over_idx)+1, 2 );

idxs(1, 1) = 1;
idxs( 2 : end, 1 )   = data_over_idx + 50;
idxs( 1 : end-1, 2 ) = data_over_idx;

% pts now contains the start of every spike

% write into matlab format
file = fopen( 'NoSpikeIdxs.txt', 'w' );
fprintf( file, '[ ');
for i = 1 : length( data_over_idx)
    current_line = sprintf( '%d : %d,...\n', idxs( i, :) );
    fprintf( file, current_line );
end
last_line = sprintf( '%d : end_idx ];', idxs( end, 1) );
fprintf( file, last_line );
fclose( file );

end