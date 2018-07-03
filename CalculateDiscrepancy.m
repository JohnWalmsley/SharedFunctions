function [ discrepancy, exp_data ] = CalculateDiscrepancy( cell, protocol, I )
%CALCULATEDISCREPANCY Summary of this function goes here
%   Detailed explanation goes here

cd ../ExperimentalData
cd(cell)
exp_data = importdata( [ protocol '_',cell,'_dofetilide_subtracted_leak_subtracted.mat' ] );
cd ..
cd ..
cd Code/
discrepancy = exp_data-I;

end

