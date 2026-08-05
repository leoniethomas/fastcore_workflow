%% Main script n°3: Models Comparison 
% This script is showing how to compare models stored inside a project.
% Models can be compared on structural, functional and sampling aspects.
% For both functional and sampling comparison, functional tests (as FBA, FBA or sampling)
% must have been previously performed (see tutorial script n°2 for single model analysis) 
% and are required to be store inside the project in the proper format.
% 
%% INITIALIZING THE ENVIRONNEMENT
initCobraToolbox();
changeCobraSolver('gurobi');
feature astheightlimit 2000;

%% LOADING PROJECT AND PARAMETERS FOR ANALYSIS
load('BRCAProject.mat');
load('defaultParametersAnalysis.csv'); % can be modified according to the user

%% COMPUTING THE COMPARISON
modelsToCompare = {"Control", "StageIV", "StageI"};
[BRCAProject, analysisIDs] = chooseActiveAnalysis(BRCAProject, modelsToCompare); % by default, the most recent analyses will be chosen 

% analysisToChoose = {'analysis_20260717_1704', 'analysis_20260717_1757', 'analysis_20260717_1841'};
% [BRCAProject, analysisIDs] = chooseActiveAnalysis(BRCAProject, modelsToCompare, analysisToChoose);

%% VISUALIZING AUTOMATICALLY GENERATED FIGURES

referenceModel = "consistentMediumConstrainedModel"; 
comparisonList = ["structuralComparison",...
    "functionalComparison", ...
    "samplingComparison"];
compID = "tutorialComparison_20260728"; 

%[BRCAProjectBis, comparisonName] = modelsComparison(BRCAProjectBis, modelsToCompare, referenceModel, compID); % only structural comp
[BRCAProject, comparisonName] = modelsComparison(BRCAProject, modelsToCompare, referenceModel, compID, comparisonList);

%% GENERATING FIGURES ON DEMAND



