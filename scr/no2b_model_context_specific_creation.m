%%

% construct context specific model based on gene expression data 

%% set up session & define input parameter file

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

% read in the parameters needed for the analysis
def_run_file = "\\atlas.uni.lux\fstc_sysbio\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model\data\def_run_paramters.txt";


%% read in all the script parameters and set working directory, directory the discretization is saved into

input_paramters = readtable(def_run_file, 'Delimiter','\t');
scr_para = cell2struct(input_paramters{:,"value"}, input_paramters{:,"slot_name"});

% add the working path to the path & set the github repo location to be the working directory
addpath(genpath(scr_para.set_working_directory));
cd(scr_para.set_working_directory);

% define a unique identifier to be used to name the different  discretizations created
date = char(datetime('now', 'Format', 'yyyyMMdd_hhss')); % to name the model and all the output 

% creating a directory for the following discretzation to be stored into - copying the defined parameters 
mkdir((scr_para.save_disc_data_to), date)
mkdir((scr_para.QC_figures_path), date)
copyfile(def_run_file, ...
         string(scr_para.save_disc_data_to) + date+ filesep + date + "_def_run_paramters.txt")
     

clear def_run_file input_paramters 



%% load preprocessed gene expression data

disc_data = "\\atlas.uni.lux\FSTC_SYSBIO\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model\discretization\20250716_1109\20250716_1109_disc_data.mat";

load(disc_data)


%% create medium for the models

med = medium(["RPMI1640.tsv", "FBS.tsv"],"./data/media");

med = med.read_medium_files(["Concentration_M","Concentration_M"])
med.medium_composition.Flux_mmol_gDW_h = - conc2Rate(med.medium_composition.Concentration_M,1e5,24,400e-12);

needed_mets = ["o2[e]", "co2[e]", "h2o[e]","h[e]", "oh1[e]"];
med = med.add_additional_rxns_boundaries(string(split(scr_para.unwanted_uptakes_export_ub, ";"))',...
                                         string(split(scr_para.unwanted_uptakes_export_lb, ";"))',...
                                         needed_mets,[])

%% start a fastcore experiment & test the consistency of input model

load(scr_para.model_used)
load(scr_para.gene_dic_file)
model_orig = model;

% BUILD generic CONSISTENT model - fast consistency check (fastcc)
exp = fastcore_experiment(model_orig,dico)
                                     
%% BUILD medium-constrained CONSISTENT model - fast consistency check (fastcc)

exp = exp.medium_constrain(med,"set_fluxes",1)


%% BUILT transcriptomics constrained CONTEXT SPECIFIC MODELS -> reconstruction using rFASTCORMICS

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


