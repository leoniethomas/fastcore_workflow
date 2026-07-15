%% Main script n°3: Models Comparison 
% This script is showing how to compare models stored inside a project.
% Models can be compared on structural, functional and sampling aspects.
% For both functional and sampling comparison, functional tests must have
% been previously performed and are required to be store inside the project
% in the proper format.
% 
%% INITIALIZING THE ENVIRONNEMENT
initCobraToolbox();
changeCobraSolver('gurobi');
feature astheightlimit 2000;

%% LOADING PROJECT AND PARAMETERS FOR ANALYSIS
load('BRCAProject.mat');
load('defaultParametersAnalysis.csv'); % can be modified according to the user

%% COMPUTING THE COMPARISON
modelsToCompare = ["Control", "StageI", "StageIV"];
[milano_project, analysisID] = chooseActiveAnalysisForComparison(milano_project, modelList, 1); % use of loopless


%% VISUALIZING AUTOMATICALLY GENERATED FIGURES

%% GENERATING FIGURES ON DEMAND