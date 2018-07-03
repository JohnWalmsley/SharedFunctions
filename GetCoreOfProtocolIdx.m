function idx = GetCoreOfProtocolIdx( protocol )
%GETCOREOFPROTOCOLIDX Summary of this function goes here
%   Detailed explanation goes here

if strcmp( protocol, 'sine_wave' )
    idx = 30000 : 65000;
elseif strcmp( protocol, 'ap' )
    idx = 5700 : 73240;
elseif strcmp( protocol, 'original_sine' )
    idx = 30000 : 80000;
elseif strcmp( protocol, 'equal_proportions' )
    idx = 30000 : 80000;
else
    disp( 'Protocol not in GetCoreOfProtocolIdx.m, Returning sine wave idxs' )
    idx = 30000 : 65000;
end

end

