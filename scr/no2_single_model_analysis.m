%% Main script n°2: Single Model Analysis 
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
% Finally, an analysis can be integrated into an already existing analysis
% thanks to addAnalysisToExistingOne.m function.
% 
%% INITIALIZING THE ENVIRONNEMENT
initCobraToolbox();
changeCobraSolver('gurobi');
feature astheightlimit 2000;

%% LOADING PROJECT AND PARAMETERS FOR ANALYSIS
load('BRCAProject_kldresults18082026.mat');
defaultParametersAnalysis = readInParamTable('defaultParametersAnalysis.csv');

%% PERFORMING ANALYSIS
% can take time depending on what's asked (especially (loopless) sampling))
wantedAnalyses = {'kld','sampling', 'FBA', 'FVA'};
analyzedModels = {'Control', 'StageI'};
BRCAProject = singleModelAnalysis(BRCAProject, defaultParametersAnalysis, analyzedModels, wantedAnalyses,1,1);

save 20260825_allSingleAnalysisLT.mat BRCAProject

%% GENERATING A REPORT
% List of wanted pathways for the report
pathwaysOfInterest = {'Citric acid cycle', 'Pyruvate metabolism', 'Glutamate metabolism', 'Alanine and aspartate metabolism'};
% List of metabolites
metsOfInterest = {'glc_D', 'pyr', 'lac_L', 'lac_D', 'gln_L', 'glu_L', 'ala_L'};

writeAnalysisReport(BRCAProject, 'StageIV', 'analysis_20260727_1852', ...
    'pathwaysOfInterest', pathwaysOfInterest, 'metsOfInterest', metsOfInterest);

%% ADDING AN ANALYSIS TO AN EXISTING ONE
BRCAProject = addAnalysisToExistingOne(BRCAProject, defaultParametersAnalysis, 'Control', 'kld', 'analysis_20260812_2130');



