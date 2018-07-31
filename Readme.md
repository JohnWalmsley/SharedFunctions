# Shared functions for discrepancy analysis

This repository contains files common to different projects in the Discrepancy Analysis projects to avoid duplication. Many of these files have been adapted from the [sine-wave repository](https://github.com/mirams/sine-wave) provided by [Beattie et al](https://physoc.onlinelibrary.wiley.com/doi/abs/10.1113/JP275733). To use the code in the projects [DiscrepancyCorrelations](https://github.com/JohnWalmsley/DiscrepancyCorrelations) and [FittingAlgorithm](https://github.com/JohnWalmsley/FittingAlgorithm), you will need to add this folder and its subfolder to the Matlab path.

## Models


#### Mex files

This folder contains Mex code for the following models, to be compiled with CVODE as follows:

[MexHH.c](Models/MexHH.c): The Hodgkin-Huxley model as used by Beattie et al

[MexMH.c](Models/MexMH.c): The Hodgkin-Huxley model in mh formulation. Note that this model will produce the same results as MexHH.c in the mh format, but is used to produce m^2h etc models as well.

[MexWang.c](Models/MexWang.c): The Wang model as used in the data supplement of Beattie et al.

There are also two Mex files that are used that additionally calculate parameter sensitivities. These use the CVODE package CVODES to compute sensitivities quickly.

[MexHHSens.c](Models/MexHHSens.c): The Hodgkin-Huxley model as used by Beattie et al.

[MexMHSens.c](Models/MexMHSens.c): The Hodgkin-Huxley model in mh formulation.

Compiled versions of all mex files are also included. If any changes to Mex files are made these must be recompiled using:

`mex -I/path_to_CVODE/include -L/path_to_CVODE/lib -lsundials_nvecserial -lsundials_cvode MexFilename.c`

or

`mex -I/path_to_CVODE/include -L/path_to_CVODE/lib -lsundials_nvecserial -lsundials_cvodes MexFilenameSens.c`

as appropriate.


#### Checking Solutions

[TestAnalyticSolutions](Models/TestAnalyticSolutions.m): allows comparison of the sensitivities and model variables from simulation and theory from the HH model. Note that the root of [SharedFunctions](https://github.com/JohnWalmsley/SharedFunctions) must be added to the path for this function to run. The remaining .m and .mat files in [Models](Models/) are called by [TestAnalyticSolutions](Models/TestAnalyticSolutions.m) only.


## Main Folder

#### Processing data

[CalculateDiscrepancy.m](CalculateDiscrepancy.m): Given string identifying a cell and a protocol, and a vector containing a (simulated) current, calculates the discrepancy between simulation and experiment.

[CalculateImposedVoltage.m](CalculateImposedVoltage.m): Given a string identifying a protocol and a time constant, generates a modified protocol that represents the delay in the amplifier and so in theory should allow a more accurate comparison between simulation and experiment. In practice these effects were negligible and so were not pursued as a source of discrepancy (although they can  be important for INa!).

[CorrectLeakCurrent.m](CorrectLeakCurrent.m): Removes the leak current from an experimental data trace. This file is no longer used as this approach can lead to strange behaviour in the activation step. See [LeakSubtractor.m](https://github.com/JohnWalmsley/Kylie-Sine-Wave-Data/blob/master/LeakSubtractor.m) for a reproduction of Kylie Beattie's approach.

[DetermineNoSpikeIdxs.m](DetermineNoSpikeIdxs.m): Helper function to find the indices close to the capacitative spikes in data. Used to write [GetNoSpikeIdx.m](GetNoSpikeIdx.m).

#### Get Experiment and Simulation Properties

[GetCoreOfProtocolIdx.m](GetCoreOfProtocolIdx.m): Given a string identifying a protocol (e.g. `'sine_wave'`), returns the indices of that protocol corresponding to its 'core', for example only the sine wave generated part of the Sine Wave protocol, without the leak or activation / deactivation steps at the beginning and end.

[GetNames.m](GetNames.m): For a cell of predictors as used in [PredictDiscrepancyUsingLinearModel.m](https://github.com/JohnWalmsley/DiscrepancyCorrelations/blob/master/Code/PredictDiscrepancyUsingLinearModel.m), this file will return the names of all predictors in the config cell.

[GetNoSpikeIdx.m](GetNoSpikeIdx.m): For a given protocol, returns all indices that are NOT near a capacitative spike. In practice this is everything except an interval of 5ms around a sudden change in potential.

[GetTemperature.m](GetTemperature.m): For a given cell identifier (e.g. `'16713110'`), returns the temperature (in **Celcius**) at which that experiment was performed.

[modeldata.m](modeldata.m): Given a string identifying a model (e.g. `'hh'`) returns default parameters (not fitted to our experimental data) and an integer identifying the model for when it is run.

#### Running Simulations

[CalculateJacobianMatrix.m](CalculatingJacobianMatrix.m): Given a voltage, model parameters, and a model type identifier from [modeldata.m](modeldata.m), calculates the Jacobian of an HH or an mh model. Used for setting sensitivity intiial conditions in [SimulatingDataSens.m](SimulatingDataSens.m).

[CalculateRates.m](CalculatingRates.m): Given a voltage or vector of voltages, model parameters, and a model type identifier from [modeldata.m](modeldata.m), calculates the rates from the model. Called by [CalculateVariables.m](CalculateVariables.m).

[CalculateSensitivityRhsVectors.m](CalculateSensitivityRhsVectors.m): Given a voltage, values of state variables, model parameters, and a model type identifier from [modeldata.m](modeldata.m), calculates the value of the RHS of the sensitivity differential equations. Used for setting sensitivity intiial conditions in [SimulatingDataSens.m](SimulatingDataSens.m).

[CalculateVariables.m](CalculateVariables.m): Given a protocol to simulate (string), a cell identifier (string), fitting protocol (cell or string) and model (string), this file uses [SimulatingDataSens.m](SimulatingDataSens.m) to perform a simulation of the protocol for the parameters determined by the fitting protocol used for the model. The function returns the following:
1. variables (Xi)
1. current (I)
1. rates (Kij)
1. fluxes (Kij Xi)
1. sensitivities ( dXi/dPj )
1. dIdp ( G dO/dp (V - Ve) )
1. voltage (V)
1. dV/dt
1. time (t)
1. conductance (G)
These variables are used for predicting discrepancy in [DiscrepancyCorrelations](https://github.com/JohnWalmsley/DiscrepancyCorrelations/).

[FindingBestFits.m](FindingBestFits.m): Given a model (string), protocol or set of protocols (cell or string), and cell identifier, returns the best fits in ascending order for that run of [CMA-ES](https://github.com/JohnWalmsley/FittingAlgorithm/blob/master/Code/cmaes.m).

[FindingBestFitsAfterMCMC.m](FindingBestFitsAfterMCMC.m): Given a model (string), protocol or set of protocols (cell or string), and cell identifier, returns the MCMC chain and likelihood values for an [MCMC](https://github.com/JohnWalmsley/FittingAlgorithm/blob/master/Code/AdaptiveMCMCStartingBestCMAES.m) run.

[SimulatingData.m](SimulatingData.m): Runs a simulation directly using the CVODE version of any model in the project. Note that, because voltage is clamped, open probabilities for models of the form m^k.h are generated from the output of [MexMH.c](Models/MexMH.c).

[SimulatingDataSens.m](SimulatingDataSens.m): Runs a simulation directly using the CVODES version of the HH or MH model in the project. Note that, because voltage is clamped, open probabilities for models of the form m^k.h are generated from the output of [MexMHSens.c](Models/MexMHSens.c). Sensitivities for the current in these models is calculated using the chain rule. Called by [CalculatingVariables.m](CalculatingVariables.m).



