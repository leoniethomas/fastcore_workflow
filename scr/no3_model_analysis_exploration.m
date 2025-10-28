%% Model analysis script 
% This script performs the analysis of a set of models.
% Overview: 
% 1. creates a model analysis object, storing all the analysis results 
% 

%% Preparations for the analysis workflow

% clear working space
clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

% add rFASTCORMICS to the path variable as well as 
addpath(genpath("C:\Users\leonie.thomas\rFASTCORMICS"))
% add scFASTCORMICS to the path -> model_analysis uses the function
% removeUnusedGenesFastbox.m -> which is a lot faster than removeUnusedGenes.m
% this is the reason why we need to load scFASTCORMICS here 
% we could also think about just adding the function to rFASTCORMICS
addpath(genpath("C:\Users\leonie.thomas\scFASTCORMICS"))
% set the solver you want to use 
changeCobraSolver("ibm_cplex");

% define script PARAMETERS

% define which of the models in your consistent model directory you want to read in 
model_id = "20251028_0709";
% do you work on a fastcore experiment object, which is created by the
% scripts in this workflow, or are you reading in your own models ? 
work_on_fastcore_exp = 1;
% path to your working directory - where your scr folder is located as well
project_path = "/Users/leonie.thomas/Documents/fastcore_workflow";
%project_path = "/Volumes/FSTC_SYSBIO/0- UserFolders/Leonie.THOMAS/projects/20250225_glynn_bulk_metabolic_model";
path_to_model_to_analyse = project_path + filesep + "context_specific_models" + filesep + model_id;
cd (project_path)
addpath(genpath(project_path))

% reading in the script parameters defined in the def_run_parameters.txt
% file - can be found in the data folder and adjusted 
input_paramters = dir(path_to_model_to_analyse + filesep + "*def_run_paramters.txt");
input_paramters = [input_paramters.folder filesep input_paramters.name];
input_paramters = readtable(input_paramters);
scr_para = cell2struct(input_paramters{:,"value"}, input_paramters{:,"slot_name"});
scr_para.model_to_load = model_id + "_cond_models.mat";
scr_para.model_workspace_to_load = model_id + "_workspace_cond_models.mat";


%% LOADING CONTEXT SPECIFIC MODELS
% There are two options. Either: 
% - you read in your own context specific models - work_on_fastcore_exp = 0
% - or you work on a fastcore_experiment object created with the previous
% workflow scripts! work_on_fastcore_exp = 1

% in case you work on your own context specific models. There are a couple
% of things to make sure, before being able to run these scripts 
% - you need to put your different models into a named struct object 
% - each model needs to contain the AA/retained_reaction vector as one
% ofthe structure slots 
% - the struct with models needs to be named 
% - you need the original model the context specific models were created on
% - and you need the dico dataframe, which is the dataframe giving the
% symbol to entrez id conversion 

if work_on_fastcore_exp
    load(path_to_model_to_analyse + filesep +   model_id + "_fastcore_exp.mat") % load the condition specific models created with rFASTCORMICS
else
    load(path_to_model_to_analyse + filesep +   model_id + "_cond_models.mat")
    load(path_to_model_to_analyse + filesep +   model_id + "_workspace_cond_models.mat")
    
    exp = fastcore_experiment(model_orig,dico);
    exp.condition_models = condition_models;
end

scr_para.results_path = project_path + filesep + "analysis" + filesep + model_id;
scr_para.objective = 'biomass';
scr_para.remove_unused_genes = 1;

mkdir(scr_para.results_path);
results = struct();

% export models to different formats
%writeCbModel(condition_models.MDA_MB231_Cont_NO,'format', 'json','fileName','model_Cont_NO.json')
%writeCbModel(condition_models.MDA_MB231_Cont_VC,'format', 'json','fileName','model_Cont_VC.json')

%% 
analysis_results = model_analysis(exp); 

%% plot jaccard similarity score for rxn presence 

[fig,J] = get_jaccard_similarity(analysis_results)

saveas(fig,scr_para.results_path + "\rxn_occurence_jaccard_distance.png");
results.jaccard = J;
clear J


%% visualize intersections rxn presence 

% get indices of outer and intersections
idx = analysis_results.get_intersection_plot("reaction_presence");

%% Pathway analysis

% get count of rxn per pathway per model -> pathway presence
analysis_results = analysis_results.get_pathway_activity(exp);
% visualize pathway presence
analysis_results.get_pathway_plot(1:15,"")

% get the subsystems for the outersection rxn

%% scatter plot based on pathway activity 

models_to_compare = string(analysis_results.model_names(1:2));
plotting_data = analysis_results.get_pathway_plot(1:200,"relative");
labeling_data = analysis_results.get_pathway_plot(1:20,"relative");

figure
scatter(plotting_data.(models_to_compare(1)), plotting_data.(models_to_compare(2)))
hold on 
text(labeling_data.(models_to_compare(1)), labeling_data.(models_to_compare(2)),...
                         labeling_data.Properties.RowNames, ...
                         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                     
%% FBA 

exp.condition_models = structfun(@(x) changeObjective(x,'biomass_reaction'),...
                             exp.condition_models,'UniformOutput',false);
exp.fba        = structfun(@(x) optimizeCbModel(x,'max','zero'),...
                                 exp.condition_models,'UniformOutput',false);
                         
fba_flux_matrix = exp.join_fba_output();

ex_rxns = exp.original_model.rxns(find(findExcRxns(exp.original_model)));
fba_exchange = fba_flux_matrix(ex_rxns(find(ismember(ex_rxns, fba_flux_matrix.Properties.RowNames))),:);
fba_exchange = fba_exchange(find(sum(abs(fba_exchange{:,:}),2) > 0),:)
fba_exchange.rxns = fba_exchange.Properties.RowNames;

writetable(fba_exchange,scr_para.results_path + filesep +  "consumption_export.xlsx") 

                         
disp("display jaccard distance based on which metabolites are ...")
disp("consumed:")
J = 1-squareform(pdist((fba_exchange{:,1:end-1} > 0)','jaccard'))
disp("and exported:")
J = 1-squareform(pdist((fba_exchange{:,1:end-1} < 0)','jaccard'))    


analysis_results.get_fba_plot(fba_flux_matrix)


%% Flux Variability Analysis

[exp.fva.minFlux, exp.fva.maxFlux] = structfun(@(x) fluxVariability(x),...
                                               exp.condition_models,'UniformOutput',false);
 
[exp.fva.minFlux,exp.fva.maxFlux] = join_fva_output(exp);

% compute FVAsimilarity

exp = exp.compute_fva_similariy()
analysis_results.get_fva_sim_plot(exp.fva_similarity)
% not sure what to do with the rxn wise fva similarity...heatplot ? not
% sure


%% run gene essentiality analysis

exp.condition_models = structfun(@(x) changeObjective(x,'biomass_reaction'),...
                             exp.condition_models,'UniformOutput',false);
                         
[grRatio, grRateKO, ...
 grRateWT, ~,...
 essential_genes_del_Rxns, ~] = structfun(@(x) singleGeneDeletion(x, 'FBA', string(x.genes), 0, 0),...
                             exp.condition_models,'UniformOutput',false);
                
analysis_results = analysis_results.add_essentiality_analysis_to_analysis_obj(exp,grRatio,grRateKO,essential_genes_del_Rxns,grRateWT)

% visualize gene essentiality
analysis_results.get_essentiallity_plots(scr_para.results_path)
% jaccard distance based on essential genes
[fig, J] = analysis_results.get_jaccard_similarity_ess_genes(0.5)
% add figure that counts essential genes per pathway


%% check for enrichment of the genes in databases 

analysis_results.enrichment = analysis_results.perform_gene_enrichment(exp,0.5);

analysis_results.visualize_enrichment()
 
%% Drug deletion
% % define a list of drugs 
% load(scr_para.gene_drug_relation_file);
% DrugList = unique(GeneDrugRelations.DrugName);
% condition_models = structfun(@(x) changeObjective(x,'biomass_reaction'),...
%                              condition_models,'UniformOutput',false);
%                          
%                          
% [grRatio, grRateKO, grRateWT] = structfun(@(x) DrugDeletion(x,'FBA',DrugList),...
%                                           condition_models,'UniformOutput',false);
%                                       
% drug_deletion_res = struct("grRatio",struct2array(grRatio),...
%                            "grRateKO",struct2array(grRateKO),...
%                            "grRateWT",struct2array(grRateWT));
% 
% drugidxs_with_an_effect = find(sum(drug_deletion_res.grRatio,2)<(size(drug_deletion_res.grRatio,2) - 0.0001));
% fig = fun.plot_clustergram(double(drug_deletion_res.grRatio(drugidxs_with_an_effect,:)),...
%                      DrugList(drugidxs_with_an_effect)',model_names,...
%                      {'Drugdeletion - grRatio KO/WT'},...
%                      [100 100 800 600],...
%                      altcolor);
% saveas(fig,scr_para.results_path + "\drug_deletion_growthRateKO_WT.png");
% 
% results.drug_deletion = drug_deletion_res;
% 
% clear DrugList grRatio grRateKO grRateWT drug_deletion_res drugidxs_with_an_effect




