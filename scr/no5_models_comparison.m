%%%% Model Comparison script - run after the Model analysis script!!!

% This script is structured in 2 parts -> structural comparison, and
% qualitative comparison
% structural analysis investigates the model structure 
% and functional analysis investigates the output of fba, fva, sampling etc

%% Set up 

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names
addpath(genpath("C:\Users\leonie.thomas\rFASTCORMICS"))
addpath(genpath("C:\Users\leonie.thomas\scFASTCORMICS"))
changeCobraSolver("gurobi");

% load the project object 
working_path = "/Users/leonie.thomas/Documents/fastcore_workflow_with_vanille";
cd (working_path)
addpath(genpath(working_path))


% read in project object
%load(working_path + filesep + "context_specific_models" + filesep + "20260119_1042" + filesep + "20260119_1042_project.mat");
load(working_path + filesep + "context_specific_models" + filesep + "20260119_1042" + filesep + "project_23012026_1453_from_vanille_fba_fva.mat");
%% First the STRUCTURAL comparison of the choosen models

% structureComparison executed on the project generating the structural
% report for the models choosen to compare 


% structureComparison function does:
%   + check if the models specified have all the needed fields ->
%     checkRequiredFieldsForModelComparison
%   + then it returns the 

%% For every analysis you need to define a list of models first which are meant to be compared 

list_model_names = ["KO","PLV", "WT"];
reference_model = "orig_model";

project = modelsComparison(project,list_model_names,reference_model,["modelStructureComparison"],"");
% modelsComparison is the function that is run once to return a set of
% predefined figures, outputs to get an overall picture of how the models
% look like in relation to each other, but the user will also be able to
% depending on this more general report do more in depth invesitgation 

%% Then Functional Comparison

% Steps of the functional comparison 
% - size of the models
% - presence in the model per subsystem 
%   - once in absoute measure -> # oversections / outersections
%   - Jaccard distance 
%   - be able to specify different subsystems in the venn diagramm
%   - then visualize all of that on a heatmap with the Jaccard distance per
%   subsystem

% show size of the models

% show compute model presence in comparison to the input model


%% Investigations up to the user (for example for specific subsystems)


% run a venn for a specific subsystem
choosen_subsystem = "Glycolysis/gluconeogenesis";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(reference_model).model.subSystems));
% create the rxns presence table from it
subsystem_feature_presence = project.comparisons.KO_vs_PLV_vs_WT.rxn_mapping_table{idx_subsystem_reference_model,:} ~= 0;
fig = plotFlexibleVenn(subsystem_feature_presence,...
                 project.comparisons.KO_vs_PLV_vs_WT.modelNames, ... 
                 "Structural model comparison: rxns presence in the " + choosen_subsystem);

choosen_subsystem = "Exchange rxns";
idx_subsystem_reference_model = find(findExcRxns(project.models.orig_model.model));
subsystem_feature_presence = project.comparisons.KO_vs_PLV_vs_WT.rxn_mapping_table{idx_subsystem_reference_model,:} ~= 0;
plotFlexibleVenn(subsystem_feature_presence,...
                 project.comparisons.KO_vs_PLV_vs_WT.modelNames, ... 
                 "Structural model comparison: rxns presence in the " + choosen_subsystem);









