%% define VARIABLES 

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names



%% define script parameters

model_id = "20250525_0950"; % check which models are available in the analysis folder 
project_path = "\\atlas.uni.lux\FSTC_SYSBIO\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model";
path_to_model_to_analyse = project_path + filesep + "context_specific_models" + filesep +model_id;
path_to_sampling_data = "." +filesep + "analysis" + filesep +model_id + filesep + "sampling" + filesep ;
cd (project_path)
addpath(genpath(project_path))
sampling_files = string(ls(path_to_sampling_data));
sampling_files = sampling_files(contains(sampling_files, "amplingResults"));
sampling_files = regexprep(sampling_files, " ", "")



%% load the created models with their whole workspace

load(path_to_model_to_analyse + "\" +   model_id + "_workspace_cond_models.mat") % load the condition specific models created with rFASTCORMICS

input_paramters = dir(path_to_model_to_analyse + "\" + "*def_run_paramters.txt");
load(scr_para.model_used) 
model_orig = model;
clear model 

input_paramters = [input_paramters.folder '\' input_paramters.name];
input_paramters = readtable(input_paramters);
scr_para = cell2struct(input_paramters{:,"value"}, input_paramters{:,"slot_name"});
scr_para.results_path = project_path + "\analysis\" + model_id;


scr_para.model_to_load = model_id + "_cond_models.mat";
scr_para.model_workspace_to_load = model_id + "_workspace_cond_models.mat";
scr_para.objective = 'biomass';
scr_para.remove_unused_genes = 1;
scr_para.gene_drug_relation_file = './data/GeneDrugRelations.mat';


altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
                       255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar
                   
%%
exp = fastcore_experiment(sampling_files,1);

exp = exp.join_fluxsum_output();
exp = exp.join_sampling_output();

%%
exp = exp.change_model_labels(["samplingResults_MDA_MB231_Cont_NO_model_20250602_090252",...
                               "samplingResults_MDA_MB231_Cont_NO_model_20250602_090543",...
                               "samplingResults_MDA_MB231_Cont_VC_model_20250602_090916",...
                               "samplingResults_MDA_MB231_Cont_VC_model_20250602_091221",...
                               "samplingResults_MDA_MB231_HERVK_C_NO_model_20250602_091549",...
                               "samplingResults_MDA_MB231_HERVK_C_NO_model_20250602_091904",...
                               "samplingResults_MDA_MB231_HERVK_C_VC_model_20250602_092222",...
                               "samplingResults_MDA_MB231_HERVK_C_VC_model_20250602_092545"],...
                              ["Cont_NO", "Cont_NO","CTR", "CTR", "HERVK_NO", "HERVK_NO", "HERVK", "HERVK"]);

%% perform PCA

exp = exp.visualize_sampling(0,... 
                             1,2,...
                             1,...
                             "samples",0);
                         
exp.visualize_sampling_rxn_distribution(3316)


%% perform differential testing - wilcoxon rank

exp.diff_flux_testing(["Cont_NO",...
                       "CTR"],1,1,scr_para.results_path + filesep)
                       

















   



