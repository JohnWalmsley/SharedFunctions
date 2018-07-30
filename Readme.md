# Shared functions for discrepancy analysis based on Sine Wave data

This folder contains files common to different projects in the Discrepancy Analysis projects to avoid duplication. Many of these files have been adapted from the sine-wave repository provided by [Beattie et al](https://physoc.onlinelibrary.wiley.com/doi/abs/10.1113/JP275733). To use the code in the projects DiscrepancyCorrelations and FittingAlgorithm, you will need to add this folder and its subfolder to the Matlab path.

## Models

This folder contains Mex code for the following models, to be compiled with CVODE as follows:

The Hodgkin-Huxley model as used by Beattie et al:
[MexHH.c](Models/MexHH.c)

The Hodgkin-Huxley model in mh formulation:
[MexMH.c](Models/MexMH.c)
(Note that this model will produce the same results as MexHH.c in the mh format, but is used to produce m^2h etc models as well)

The Wang model as used in the data supplement of Beattie et al.
[MexWang.c](Models/MexWang.c)

There are also two Mex files that are used that additionally calculate parameter sensitivities. These use the CVODE package CVODES to compute sensitivities quickly.

The Hodgkin-Huxley model as used by Beattie et al:
[MexHHSens.c](Models/MexHHSens.c)

The Hodgkin-Huxley model in mh formulation:
[MexMHSens.c](Models/MexMHSens.c)

Compiled versions of all mex files are also included. If any changes to Mex files are made these must be recomplied using:

`mex -I/path_to_CVODE/include -L/path_to_CVODE/lib -lsundials_nvecserial -lsundials_cvode MexFilename.c`

Or

`mex -I/path_to_CVODE/include -L/path_to_CVODE/lib -lsundials_nvecserial -lsundials_cvodes MexFilenameSens.c`

as appropriate.

[TestAnalyticSolutions](Models/TestAnalyticSolutions.m) allows comparison of the sensitivities and model variables from simulation and theory from the HH model. Note that the root of [SharedFunctions](/) must be added to the path for thius function to run. The remaining .m and .mat files in [Models](Models/) are called by [TestAnalyticSolutions](Models/TestAnalyticSolutions.m) only.





