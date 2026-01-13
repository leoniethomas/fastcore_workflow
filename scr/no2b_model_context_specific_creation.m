%%

% construct context specific model based on gene expression data 

%% set up session & define input parameter file

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

% read in the parameters needed for the analysis
%def_run_file = "/Volumes/FSTC_SYSBIO/0- UserFolders/Leonie.THOMAS/projects/20250225_glynn_bulk_metabolic_model/data/def_run_paramters.txt";

% define which discretization run you want to use for model building
disc_data_id = "20251125_0837";
def_run_file = "/Users/leonie.thomas/Documents/GSNOR_metab_models/fastcore_workflow/discretization" + filesep + ...
                disc_data_id + filesep + disc_data_id + "_def_run_paramters.txt";


%% read in all the script parameters and set working directory, directory the discretization is saved into

scr_para = read_in_run_def_file(def_run_file);

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
changeCobraSolver("gurobi")

%% load preprocessed gene expression data

disc_data = string(scr_para.set_working_directory) + filesep + "discretization" + filesep + disc_data_id + filesep + disc_data_id + "_disc_data.mat";

load(disc_data)


%% create medium for the models

med = medium(["DMEM.tsv", "FBS.tsv"],"./data/media");

med = med.read_medium_files(["Concentration_M","Concentration_M"]); 
med.medium_composition.Flux_mmol_gDW_h = - conc2Rate(med.medium_composition.Concentration_mM,...
                                                       scr_para.cell_conc_cells_per_ml,...
                                                       scr_para.t_in_hours,...
                                                       scr_para.cell_dry_weight_gDW); % TODO put the values from the description file here 

med = med.add_additional_rxns_boundaries(string(split(scr_para.unwanted_uptakes_export_ub, ";"))',...
                                         string(split(scr_para.unwanted_uptakes_export_lb, ";"))',...
                                         string(split(scr_para.needed_mets_medium, ";"))',...
                                         [])

%% start a fastcore experiment & test the consistency of input model

model = load(scr_para.model_used);
model = model.(string(fieldnames(model)));
model = generateRules(model);
load(scr_para.gene_dic_file)
model_orig = model;

% BUILD generic CONSISTENT model - fast consistency check (fastcc)
exp = fastcore_experiment(model_orig,dico)
                                     
%% BUILD medium-constrained CONSISTENT model - fast consistency check (fastcc)

%% TODO -> set the concentraiton as fluxes !!
exp = exp.medium_constrain(med,"set_fluxes",0)


%% BUILT transcriptomics constrained CONTEXT SPECIFIC MODELS -> reconstruction using rFASTCORMICS

exp.data = data;
exp.data.source = disc_data;

%%%%%%%%%%%%%%%%%%%%%%  set settings used for the fastcormics run 
% the unpenalized argument is currently not part of fastcormics4cobra_v2
%optional_settings.unpenalized = model_orig.rxns(ismember(vertcat(model_orig.subSystems{:}), ... 
%                                                         strsplit(scr_para.unpenalizedSystems,";")));
% forcing the medium in by setting it into fun option
biomass_rxn = {'biomass_reaction'}  % 'biomass_reaction' for recon3d

optional_settings.func = {'DM_atp_c_', biomass_rxn{:,:},med.medium_composition.ExRxns_Recon3D{:}}; %biomass_maintenance %-> c
% composition constrain the model by setting the optional setting to the medium composition
optional_settings.medium = med.medium_composition.Mets_Recon3D; %(add media instead)
optional_settings.not_medium_constrained = scr_para.not_medium_constrained;

%%%%%%%%%%%%%%%%%%%%%%  
                        
condition_models = struct();

condition_column = scr_para.columns_to_define_model_samples_on;
% get the index of the samples in every defined group
for cond = unique(data.metadata.(condition_column))'
         % transform the array, the for loop loops over the rows, so if the elements over which you want to loop over are defined in cells in one row /not column then the for loop will concat all elements instead of looping over them
        idx = contains(data.metadata.(condition_column),cond);

        disp("condition for which the samples are filtered: " + cond + newline + " ----------------####################### ------------------------");
        
        % run rfastcormics on consistent global metabolic model
        tic; % mearuse the time the model takes to run
        [model_cond,retainedRxns, indicesCompletedCoreOrig] = rFastcormics4cobra_v2(exp.medium_constrained_model,data.discretized(:,idx), ...
                                                                                    cellstr(data.feature_names_norm), dico,...
                                                                                    scr_para.consensus_proportion, scr_para.epsilon,...
                                                                                    optional_settings, biomass_rxn, 0, 0);
        model_cond.running_time = toc;
        model_cond.used_data = data.discretized(:,idx); % add the data used for the model to the resulting model
        model_cond.sample_metadata = data.metadata(idx,:); % add metadta of the samples used to compute the model!
        model_cond.retainedRxns = retainedRxns;
        model_cond.indicesCompletedCoreOrig = indicesCompletedCoreOrig;
        condition_models.(strrep(cond{:},"-","_")) = model_cond;
        if ~sum(contains(model_cond.rxns,"biomass_reaction")) == 1
            error("Biomass_rxn is lost!!")
        end
        if length(model_cond.rxns)~=length(fastcc(model_cond, scr_para.epsilon, 0, 0, 'original'))
            error("Model generated by fascormics is not consistent!")
        end
        disp("number of samples for which this condition was modelled on: " + size(data.discretized(:,idx),2) + newline + " ----------------####################### ------------------------")
end

exp.optional_settings = optional_settings;
exp.condition_models = condition_models;



%%

clear A AA idx xi x TXT tsquared model_cond condition_column cond
cd(scr_para.set_working_directory)
mkdir(fullfile(scr_para.save_models_to, date))

exp_file_name = fullfile(scr_para.save_models_to,date,[ date '_fastcore_exp.mat']);  % Convert datetime object to string
disp(exp_file_name);
md_file_name = fullfile(scr_para.save_models_to,date,[ date '_cond_models.mat']);  % Convert datetime object to string
disp(md_file_name);
dat_file_name = fullfile(scr_para.save_models_to,date,[ date  '_workspace_cond_models.mat']);  % Convert datetime object to string
disp(dat_file_name);

save(md_file_name, 'condition_models')
save(dat_file_name)
save(exp_file_name, 'exp')

