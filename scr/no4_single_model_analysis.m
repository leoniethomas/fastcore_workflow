%% Single model analysis script 
% This script is showing how to perform individual analysis on context-specific models stored in a project object.
% The first function, singleModelAnalysis.m, allows to perform on each model the following analysis:
%   - FBA
%   - FVA
%   - sampling
%   - looplessSampling (using cycleFreeFlux, from Desouki et al., 2015, COBRA version)
%   - FDR correction for sampling (https://doi.org/10.1016/j.jbi.2024.104597)
%   - singleGeneDeletion
%   - doubleGeneDeletion
% Parameter choice of each wanted analysis must be provided in a table
% format following the one given in our default table (defaultParametersAnalysis.csv).
% Both parameters and performed analyses are stored in the project object as well. A second
% function, writeAnalysisReport.m, generates a report based on the
% performed analysis, including tables with FBA/FVA/flux sum for exchangers
% or given metabolites, reactions, or pathways.
% 
%% INITIALIZING THE ENVIRONNEMENT
initCobraToolbox();
changeCobraSolver('gurobi');
feature astheightlimit 2000;
%% LOADING PROJECT AND PARAMETERS FOR ANALYSIS
load('BRCAProject.mat');
load('defaultParametersAnalysis.csv'); % can be modified according to the user

%% PERFORMING ANALYSIS
% can take time depending on what's asked (especially (loopless) sampling))
wantedAnalyses = {'FBA', 'FVA'};
analyzedModels = {'Control', 'StageI', 'StageIV'};
BRCAProject_test = singleModelAnalysis(BRCAProject, defaultParametersAnalysis, analyzedModels, wantedAnalyses);

%% ADD AN ANALYSIS TO AN EXISTING ONE
BRCAProject_test = addAnalysisToExistingOne(BRCAProject_test, defaultParametersAnalysis, 'StageIV', 'singleGeneDeletion', 'analysis_20260624_1625');

%% REPORT
% List of wanted pathways for the report
pathwaysOfInterest = {'Citric acid cycle', 'Pyruvate metabolism', 'Glutamate metabolism', 'Alanine and aspartate metabolism'};
% List of metabolites
metsOfInterest = {'glc_D', 'pyr', 'lac_L', 'lac_D', 'gln_L', 'glu_L', 'ala_L'};

writeAnalysisReport(BRCAProject_test, 'StageIV', 'analysis_20260624_1625', ...
    'pathwaysOfInterest', pathwaysOfInterest, 'metsOfInterest', metsOfInterest);