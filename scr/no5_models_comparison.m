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
% load(working_path + filesep + "context_specific_models" + filesep + "20260119_1042" + filesep + "project_23012026_1453_from_vanille_fba_fva.mat");
% old_project = project;
% load(working_path + filesep + "context_specific_models" + filesep + "20260119_1042" + filesep + "project_28012026_1508_obj_vanille_sampling.mat");
% samp_project = project;
load(working_path + filesep + "context_specific_models" + filesep + "20260119_1042" + filesep + "project_23012026_1453_28012026_1508_obj_vanille_sampling.mat")
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
identifier = "";

[project,comparison_name] = modelsComparison(project,list_model_names,reference_model,["modelStructureComparison"],identifier);
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
subsystem_feature_presence = project.comparisons.KO_vs_PLV_vs_WT__.rxn_mapping_table{idx_subsystem_reference_model,:} ~= 0;
fig = plotFlexibleVenn(subsystem_feature_presence,...
                 project.comparisons.KO_vs_PLV_vs_WT__.modelNames, ... 
                 "Structural model comparison: rxns presence in the " + choosen_subsystem);

choosen_subsystem = "Exchange rxns";
idx_subsystem_reference_model = find(findExcRxns(project.models.orig_model.model));
subsystem_feature_presence = project.comparisons.KO_vs_PLV_vs_WT__.rxn_mapping_table{idx_subsystem_reference_model,:} ~= 0;
plotFlexibleVenn(subsystem_feature_presence,...
                 project.comparisons.KO_vs_PLV_vs_WT__.modelNames, ... 
                 "Structural model comparison: rxns presence in the " + choosen_subsystem);

%% genes of interest 
gene = "5831"; 
genes_in_model = string(project.models.(reference_model).model.genes);
genes_of_interest =  genes_in_model(find(contains(genes_in_model , gene)));

findRxnsFromGenes(project.models.(reference_model).model,char(genes_of_interest(1)))

%% get the rxns overview with FVA,FBA and reduced cost for a specified subsystem

choosen_subsystem = "Exchange/demand reaction";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(reference_model).model.subSystems));
get_flux_plot(project,"KO_vs_PLV_vs_WT__",idx_subsystem_reference_model, "FVA",true)



%% get rxns which show a reduced cost ~= 0 

replacement_value = "analysis.FBA.basis.reducedcost"; % get the fba solution values
ordered_reducedCost_matrix = getOrderedFeatureMatrix(project,project.comparisons.KO_vs_PLV_vs_WT__.modelNames,"rxns","orig_model",replacement_value);
reduced_cost_idx = find(sum(ordered_reducedCost_matrix,2) ~= 0);
get_flux_plot(project,"KO_vs_PLV_vs_WT__",reduced_cost_idx, ...
              "threshold_flux","all","FVA", true,'reducedCost',true ,...
              'title_plots',"Functional model comparison: all reactions with a ~= 0 reduced cost");
 
%% get rxns associated with a specific metabolite 

met_name_pattern = "^pro_L[*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(reference_model).model.mets, met_name_pattern, 'once')));
met_names = project.models.(reference_model).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(reference_model).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(reference_model).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","none","FVA", true ,'threshold_flux', "all", ...
              'title_plots',"Functional model comparison: all reactions including proline");

%% get rxns associated with a specific metabolite 

rxn_name_pattern = "^P5CR.*";
idx_rxn_matches = find(~cellfun(@isempty, regexp(project.models.(reference_model).model.rxns, rxn_name_pattern, 'once')));
rxn_names = project.models.(reference_model).model.rxns(idx_rxn_matches);

[~,rxns_id] = ismember(rxn_names, project.models.(reference_model).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_id, ...
              "threshold_flux","none","FVA", true ,'threshold_flux', "all", ...
              'title_plots',"Functional model comparison: all reactions including proline");

%% get rxns associated with a specific metabolite 

met_name_pattern = "^pyr[.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(reference_model).model.mets, met_name_pattern, 'once')));
met_names = project.models.(reference_model).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(reference_model).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(reference_model).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","none","FVA", true ,...
              'title_plots',"Functional model comparison: all reactions including pyruvate");
%% get rxns associated with a specific metabolite 

met_name_pattern = "^succ[.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(reference_model).model.mets, met_name_pattern, 'once')));
met_names = project.models.(reference_model).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(reference_model).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(reference_model).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","none","FVA", true ,...
              'title_plots',"Functional model comparison: all reactions including succinate");

%% get rxns associated with a specific metabolite 

met_name_pattern = "^lac.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(reference_model).model.mets, met_name_pattern, 'once')));
met_names = project.models.(reference_model).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(reference_model).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(reference_model).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","none","FVA", true ,...
              'title_plots',"Functional model comparison: all reactions including lactate");

%% SAMPLING ANALYSIS


project = modelsComparisonSampling(project,comparison_name);


%%

rxn_color = ["EX_glc_D[e]", "EX_gln_L[e]", "EX_pyr[e]","EX_so3[e]", "EX_retfa[e]"];

visualize_sampling_landscape(project,comparison_name,[1 2],rxn_color)

%%
met_name_pattern = "^lac_D[*";
idx_mets_matches = find(~cellfun(@isempty, regexp(project.models.(reference_model).model.mets, met_name_pattern, 'once')));
met_color = string(project.models.(reference_model).model.mets(idx_mets_matches));

visualize_sampling_landscape(project,comparison_name,[1 2],met_color', "fluxsum")

%% grouped boxplot for sampling results


trial1 = rand(5,7);
trial2 = rand(10,7);
trial3 = rand(15,7);

% These grouping matrices label the columns:
grp1 = repmat(1:7,size(trial1,1),1);
grp2 = repmat(1:7,size(trial2,1),1);
grp3 = repmat(1:7,size(trial3,1),1);

% These color matrices label the matrix id:
clr1 = repmat(1,size(trial1));
clr2 = repmat(2,size(trial2));
clr3 = repmat(3,size(trial3));

% Combine the above matrices into one for x, y, and c:
x = [grp1;grp2;grp3];
y = [trial1;trial2;trial3];
c = [clr1;clr2;clr3];

% Convert those matrices to vectors:
x = x(:);
y = y(:);
c = c(:);

% Multiply x by 2 so that they're spread out:
x = x*2;

% Make the boxchart, 
boxchart(x(:),y(:),'GroupByColor',c(:))

% Set the x ticks and labels, and add a legend
xticks(2:2:14);
xticklabels(1:7)
xlabel('Category')
legend(["Trial 1" "Trial 2" "Trial 3"],'Location','NorthOutside')



%% 

met_name_pattern = "^pyr[*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(reference_model).model.mets, met_name_pattern, 'once')));
met_names = project.models.(reference_model).model.mets(idx_met_matches);

[~,met_id] = ismember(met_names, project.models.(reference_model).model.mets);


% matrix with the fluxsumvalues of interest 

data = project.comparisons.(comparison_name).ordered_samples_fluxsum(met_id,:);
data = data(find(any(data ~= 0,2)),:);

samples_matrix = project.comparisons.(comparison_name).sample_model_labels;


% Your data
% data: 4 x 6000
% met_names: 4x1 cell
% samples_matrix: 1 x 6000, values like "KO","PLV","WT"

[nMet, nSamples] = size(data);

% Step 1: x-axis grouping (metabolites)
x = repmat(1:nMet, nSamples, 1); % size: nSamples x nMet

% Step 2: y-values
y = data'; % size: nSamples x nMet
y = y(:);

% Step 3: color grouping based on sample labels
% Convert sample labels to numeric codes
[~, c_numeric] = ismember(samples_matrix, unique(samples_matrix)); 
% samples_matrix is 1x6000, unique(samples_matrix) = ["KO","PLV","WT"]
% c_numeric will be 1 for KO, 2 for PLV, 3 for WT

% Repeat c_numeric for each metabolite
c = repmat(c_numeric', 1, nMet); % size: nSamples x nMet
c = c(:);

% Step 4: flatten x
x = x(:) * 2; % optional: spread boxes

% Step 5: plot
boxchart(x, y, 'GroupByColor', c);

% Step 6: x-tick labels for metabolites
xticks((1:nMet)*2);
xticklabels(met_names);

ylabel('Value');
title('Grouped Boxchart by Metabolite and Sample Type');

% Optional: add legend
legend(unique(samples_matrix));



