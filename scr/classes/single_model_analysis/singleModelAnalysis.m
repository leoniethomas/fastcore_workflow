function proj = singleModelAnalysis(proj, modelList, analyses)
% this function runs the analysis on a model or a list or models and stores
% the results in a structure.
% Inputs: Project, Names of the models, List of wanted analysis (all by
% default)
% Available analysis:
%   FBA
%   FVA
%   flux_sum
%   sampling
%   single_gene_deletion
%   double_gene_deletion
%   enrichment ?
% Output : project with a field analysis

%% Check the structure of project
% functions
% -> CheckingStructureForAnalysis(proj, modelList)
% Inside a model field, at least a model, but depends on what's necessary
% for the asked analyses
% -> checkRequiredFieldForAnalyses
%% Run analysis
% loop through the models

