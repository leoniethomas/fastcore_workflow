%%

% construct context specific model based on gene expression data 

%% Setup and parameters


clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

% read in the parameters needed for the analysis
vis_only_met_genes = 1;
def_run_file = "\\atlas.uni.lux\FSTC_SYSBIO\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model\discretization\20250525_0912\20250525_0912_def_run_paramters.txt";
disc_data = "\\atlas.uni.lux\FSTC_SYSBIO\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model\discretization\20250525_0912\20250525_0912_disc_data.mat";

input_paramters = readtable(def_run_file, 'Delimiter','\t');

% Define your array of field names
scr_para = cell2struct(input_paramters{:,"value"}, input_paramters{:,"slot_name"});

addpath(genpath(scr_para.set_working_directory));
cd(scr_para.set_working_directory);
date = char(datetime('now', 'Format', 'yyyyMMdd_hhss')); % to name the model and all the output 
mkdir((scr_para.save_models_to), date)
mkdir((scr_para.QC_figures_path), date)

copyfile(def_run_file, ...
         string(scr_para.save_models_to) + date+ filesep + date + "_def_run_paramters.txt")
     

clear def_run_file input_paramters 


%% load model and dico

load(scr_para.model_used) 
load(scr_para.gene_dic_file)
load(disc_data)

%% Load medium
% + set lower boundaries to values in input file
% + set the rest of the exchange reactions to 0
% I should still fix this part -> cause this is only applying for my input
% files 
% -> make a function out of it 

[NUM,TXT,RAW]=xlsread(scr_para.medium_used_file);
met_name = RAW(2:end,find(matches(RAW(1,:),scr_para.medium_used_naming_colum)));
flux = cell2mat(RAW(2:end,find(matches(RAW(1,:),scr_para.medium_used_concentration_column))));
rxn=cellfun(@(x)['EX_' x],met_name,'uni',false);
model.medium = table(met_name,rxn, flux);
[~,idx, idx_fluxes_in_model] = intersect(model.medium.rxn,model.rxns);   
model.lb(idx_fluxes_in_model) = -model.medium{idx,"flux"};

[EX, UPT] = findExcRxns(model);
needed_mets = ["o2[e]", "co2[e]", "h2o[e]","h[e]", "oh1[e]"];
Ex_to_close = setdiff(model.rxns(findExcRxns(model)),...
                                 [model.medium.rxn; findRxnsFromMets(model, needed_mets)]);

% set the rxns lower and upper bound, rxns that we set that we do not want
% to have, are those also reasonable in my case, for my data ? 
model.ub(find(ismember(model.rxns,split(scr_para.unwanted_uptakes_export_ub, ";"))))=0; 
model.lb(find(ismember(model.rxns,split(scr_para.unwanted_uptakes_export_lb, ";"))))=0; 

model.lb(findRxnIDs(model, Ex_to_close))=0; 

clear idx_fluxes_in_model idx
%% BUILD generic CONSISTENT model - fast consistency check (fastcc)

A = fastcc_4_rfastcormics(model, 1e-4, 1);

% remove non consistent reactions from model
model=removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));
% check if the biomass reactions are still there
model.rxns(find(contains(model.rxns,'biomass')))

model_orig=model;
% check if the created model is now really consistent
A = fastcc_4_rfastcormics(model, 1e-4, 1);

clear A model

%% BUILT CONTEXT SPECIFIC MODELS -> reconstruction using rFASTCORMICS

load(scr_para.gene_dic_file)

subSys=vertcat(model_orig.subSystems{:});
 
optional_settings.unpenalized = model_orig.rxns(ismember(subSys, ...
                                                         strsplit(scr_para.unpenalizedSystems,";")));

% you need to input the exchange rxns names, to force the model exchange
% rxns in, after check f
optional_settings.func = {'DM_atp_c_', 'biomass_reaction',model_orig.medium.rxn{:}}; %biomass_maintenance %-> c

% these two things need to be set like this otherwise the model will not be
% medium constrained!
optional_settings.medium = model_orig.medium.met_name; %(add media instead)
optional_settings.not_medium_constrained = scr_para.not_medium_constrained;

biomass_rxn = {'biomass_reaction'} 

                        
condition_models = struct();

condition_column = scr_para.columns_to_define_model_samples_on;
% get the index of the samples in every defined group
for cond = unique(data.metadata.(condition_column))'
         % transform the array, the for loop loops over the rows, so if the elements over which you want to loop over are defined in cells in one row /not column then the for loop will concat all elements instead of looping over them
        idx = contains(data.metadata.(condition_column),cond);

        disp("condition for which the samples are filtered: " + cond + newline + " ----------------####################### ------------------------");
        
        % run rfastcormics on consistent global metabolic model
        tic; % mearuse the time the model takes to run
        [model_cond,AA] = fastcormics_RNAseq(model_orig,data.discretized(:,idx), ...
                                             data.feature_names_norm, dico, biomass_rxn, str2double(scr_para.already_mapped_tag),...
                                                str2double(scr_para.consensus_proportion), str2double(scr_para.epsilon), optional_settings);
        model_cond.running_time = toc;
        model_cond.used_data = data.discretized(:,idx); % add the data used for the model to the resulting model
        model_cond.sample_metadata = data.metadata(idx,:); % add metadta of the samples used to compute the model!
        model_cond.AA = AA;
        condition_models.(strrep(cond{:},"-","_")) = model_cond;
        disp("number of samples for which this condition was modelled on: " + size(data.discretized(:,idx),2) + newline + " ----------------####################### ------------------------")
end


%%

clear A AA idx xi x TXT tsquared model_cond condition_column cond

md_file_name = [scr_para.save_models_to  date '/' date '_cond_models.mat'];  % Convert datetime object to string
disp(md_file_name);
dat_file_name = [scr_para.save_models_to   date '/' date  '_workspace_cond_models.mat'];  % Convert datetime object to string
disp(dat_file_name);

save(md_file_name, 'condition_models')
save(dat_file_name)


