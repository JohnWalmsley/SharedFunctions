function [P,V] = FindingBestFitsAfterMCMC(model,protocols,exp_ref)
% Function imports MCMC chains and likelihood values for relevant model/protocol/exp_ref combination
cd  ../MCMCResults
pathToFolder = '.';
files = dir(fullfile(pathToFolder,'*'));

P=[];
V=[];

if isa( protocols, 'cell' )
    fitting_protocol = protocols{ 1 };
    for pr = 2 : length( protocols )
        fitting_protocol = [ fitting_protocol '_' protocols{ pr } ];
    end
else
    fitting_protocol = protocols;
end

for i = 1:numel(files)
    
    k = regexp(files(i).name,['MCMCLikelihood_',exp_ref,'_',model,'_',fitting_protocol ]);

    if k>=1
        V=[V;importdata(files(i).name)];
        S=files(i).name;
        
        mS = strrep(S, 'MCMCLikelihood', 'MCMCChain');
        q=exist('P');
        if q==1
            P=[P;importdata(mS)];
            
        else
            P=importdata(mS);   
        end
    end
end

cd ..

cd Code
