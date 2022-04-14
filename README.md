## Introduction
This directory contains the codes to perform all the analyses and reproduce all the figures and tables presented in the paper: Ng HM, Jiang BY, Wong KY. Penalized estimation of a class of single-index varying-coefficient models for integrative genomic analysis. 2022.

Before running the codes, install the R-package **psivcm** which can be found in the home directory. The programs also require the R-packages  **grpreg**, **splines2**, **survival**, and **SurvC1**.

## Simulation Settings
The simulation settings considered in the paper are as follows:

* Setting 1 (`model = 1`): the simulation study in Section 3.

* Setting 2 (`model = 2`): the additional simulation study in Section S3.1 concerning a main effect model.

* Setting 3 (`model = 3`): the additional simulation study in Section S3.2 concerning a model with non-linear main effects.

* Setting 4 (`model = 4`): the additional simulation study in Section S3.3 concerning a single-index varying coefficient model where the continuous outcome variable has larger noise than Setting 1.


## Simulation Data Sets
Users can simulate the random data by running the program `SimulateData.R` in the directory `Simulation`. It will generate the initial values of single-index parameter and 101 simulation data sets (100 training data sets and a validation data set) for each simulation setting to the directory `Simulation/SimulationData`.

The simulation data sets used in the paper are available, the file `Simulation-beta0.csv` and the zip files (`SimulationData-p20.zip`, `SimulationData-p50.zip`, `SimulationData-p100.zip`, `SimulationData-p300-0-50.zip`, `SimulationData-p300-51-100.zip`, and `AdditionalSimulationData.zip`) can be found in the directory `Simulation/SimulationData`. To reproduce the analyses, download the zip files and copy the files to the directory `Simulation/SimulationData`.

The simulation data consist of the file `Simulation-beta0.csv`, 404 files with names in the form of `SimulationData-p[number of covariates in X]-[replication number].csv`, and 1212 files with names in the form of `AdditionalSimulationData-setting[setting number]-[number of covariates in X]-[replication number].csv`. For details, see the following data documentation:

* The file `Simulation-beta0.csv` stores the initial values of single-index parameter that are used for all simulation replicates. Each row in a file contains a set of initial values.

* Each of the 404 files contains the simulation data for 500 subjects (5000 subjects for `replication number = 0`) for a specific simulation replicate under setting 1. Each row in a file contains data for a subject. The first element on each row is the continuous outcome, the second element and the third element are the right-censored outcome corresponding to the observed time and the event indicator (1 for observed event time and 0 for right-censored time), respectively. The remaining elements are the observed values of X and U. 

* Each of the 1212 files contains the responses for other simulation settings (settings 2 to 4). Each row in a file contains data for a subject. The first element on each row is the continuous outcome, the second element and the third element are the right-censored outcome corresponding to the observed time and the event indicator, respectively.


## Simulation Studies

Run the program `SimulationAnalysis.R` in the directory `Simulation`. The program performs analysis on 100 simulation data sets and validates on a validation data set (with `replication number = 0`). For each simulation setting and estimation method (and initial values of single-index parameter for the proposed methods), the program generates 2 output files in the directory `Simulation/SimulationResults`. Each row of the output file contains the simulation results for a replication. Files with prefix `SimulationResults-` in the name contain a summary of the model performance, including the sensitivity, FDR, cardinality, MSE, C-index, and lambda value(s). Files with prefix `Estimates-` in the name contain the estimated parameters of the model with smallest modified BIC value.

To reproduce all the analyses in the paper, users should follow the comments stated in lines 12-19 of the program `SimulationAnalysis.R` and change the simulation setting accordingly. Since some of the programs may take long time to run, we have included all the simulation results in `SimulationResults.zip` in the directory `Simulation/SimulationResults`. Users may copy the intermediate results to `Simulation/SimulationResults`. The program `SummerizeResults.R` in the directory `Simulation` summarizes the simulation results for each setting over 100 replicates.


## Analysis of TCGA data

We provided two real data examples for the application of the proposed method. The raw data downloaded from [UCSC Xena data hubs](https://xena.ucsc.edu) can be found in the directory `RealDataAnalysis/RealDataAnalysisData`. The zip file `RealDataAnalysis-beta0.zip` contains the files `RealData-NSCLC-beta0.csv` and `RealData-LGG-beta0.csv` that store 50 initial values of single-index parameter for each analysis. Each row in a file contains a set of initial values. (Users can also generate initial values of the single-index parameter by switching `FALSE` in lines 16 and 154 of the program `RealDataAnalysis.R` into `TRUE`. The program will generate the files `RealData-NSCLC-beta0.csv` and `RealData-LGG-beta0.csv` that store 50 initial values of single-index parameter to the directory `RealDataAnalysis/RealDataAnalysisData`.)  The zip file `RealDataAnalysis-data.zip` contains all data files required for the real data analyses.

Download `RealDataAnalysis-beta0.zip` and `RealDataAnalysis-data.zip`, and copy the files to `RealDataAnalysis/RealDataAnalysisData`. Run the program `DataProcessing.R` to extract relevant data for the analyses. The processed data will be written in the same directory. Run the program `RealDataAnalysis.R` in the directory `RealDataAnalysis`. This program performs analyses on the NSCLC and LGG data sets. For each analysis, the program generates an output file in the directory `RealDataAnalysis/RealDataAnalysisResults`. For the proposed methods, each row in the output file contains the estimated regression parameters under each initial value of the single-index parameter.


## Generation of Figures

Upon completion of all analyses, run the programs `plotFigures1andS1toS7.R` in the directory `Simulation` and `plotFigures2to4.R` in the directory `RealDataAnalysis`, which generate the figures for the estimated coefficient. These programs generate Figures 1 and S1 to S7 to the directory `Simulation/SimulationResults` and Figures 2 to 4 to the directory `RealDataAnalysis/RealDataAnalysisResults`.
