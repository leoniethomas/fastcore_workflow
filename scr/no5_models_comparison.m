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
project = load(working_path + filesep + "context_specific_models" + filesep + "20260119_1042" + filesep + "20260119_1042_project.mat");

%% First the STRUCTURAL comparison of the choosen models

% structureComparison executed on the project generating the structural
% report for the models choosen to compare 


% structureComparison function does:
%   + check if the models specified have all the needed fields ->
%     checkRequiredFieldsForModelComparison
%   + 



%% Then Functional Comparison


