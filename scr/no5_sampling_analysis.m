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
sampling_files = sampling_files(contains(sampling_files, "samplingResults_MDA_MB231_Cont_"));
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

%% check fluxsum 

[~,income_flux] = exp.fastcore_runs.samplingResults_MDA_MB231_Cont_NO_model_20250602_090252.compute_flux_sum(1);
[~,outgo_flux] = exp.fastcore_runs.samplingResults_MDA_MB231_Cont_NO_model_20250602_090252.compute_flux_sum(1,0);
outgo_flux = abs(outgo_flux);

%%
exp = exp.change_model_labels(["samplingResults_MDA_MB231_Cont_NO_model_20250602_090252",...
                               "samplingResults_MDA_MB231_Cont_NO_model_20250602_090543",...
                               "samplingResults_MDA_MB231_Cont_VC_model_20250602_090916",...
                               "samplingResults_MDA_MB231_Cont_VC_model_20250602_091221"],...
                              ["Cont\_NO", "Cont\_NO","CTR", "CTR"]);

%% perform PCA

exp = exp.visualize_sampling(0,... 
                             1,2,...
                             1,...
                             "fluxsum",0);
                         
      %%                 
exp.visualize_sampling_rxn_distribution(find(matches(exp.rxn_names,"PDHm")))
exp.visualize_sampling_rxn_distribution(find(matches(exp.rxn_names,"PGK")))
exp.visualize_sampling_rxn_distribution(find(matches(exp.rxn_names,"PFK")))
exp.visualize_sampling_rxn_distribution(find(matches(exp.rxn_names,"FBA")))
exp.visualize_sampling_rxn_distribution(find(matches(exp.rxn_names,"HEX1")))
exp.visualize_sampling_rxn_distribution(find(matches(exp.rxn_names,"G6PDH2r")))


%%


exp.visualize_sampling_rxn_distribution(find(matches(exp.rxn_names,"ACONTm")),"fluxsum")
exp.visualize_sampling_rxn_distribution(find(matches(exp.rxn_names,"ACONTm")))

%%
% writeCbModel(condition_models.MDA_MB231_Cont_VC,'format','json','fileName','VC_model.json')
% writeCbModel(condition_models.MDA_MB231_Cont_NO,'format','json','fileName','NO_model.json')

%% create escher input csv 

writetable(table(exp.rxn_names, mean(exp.samples(:,1:4000),2)),"sampling_data_NO.csv")
writetable(table(exp.rxn_names, mean(exp.samples(:,4000:8000),2)),"sampling_data_VC.csv")


%% fluxes single samples

i = 123
writetable(table(exp.rxn_names, mean(exp.samples(:,i),2)),"sampling_data_ev_"+ string(i) + ".csv")
%% create escher input csv fluxsum

met_names_escher = regexprep(exp.met_names, '\[(.)\]$', '_$1');

writetable(table(met_names_escher, mean(exp.fluxsum(:,1:4000),2)),"sampling_data_fluxsum_NO.csv")
writetable(table(met_names_escher, mean(exp.fluxsum(:,4000:8000),2)),"sampling_data_fluxsum_VC.csv")
%% save  model and sample mean

writeCbModel(exp.fastcore_runs.SamplingResults_medium_1500_model_1.model,'format','json','fileName','ev_medium_1_model.json')
writetable(table(exp.fastcore_runs.SamplingResults_medium_1500_model_1.model.rxns,...
                 mean(exp.fastcore_runs.SamplingResults_medium_1500_model_1.sampling,2)),"sampling_data_medium_ev_1.csv")

%% find samples which are different to visualize them 
dist_samples = pdist(exp.samples(:,1:4000));

%% perform differential testing - wilcoxon rank

exp.diff_flux_testing(["Cont\_NO",...
                       "CTR"],1,1,scr_para.results_path + filesep)
                       

















   



