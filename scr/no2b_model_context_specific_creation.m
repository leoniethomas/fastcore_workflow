%%

% construct context specific model based on gene expression data 

%% set up session & define input parameter file

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

% define which discretization run you want to use for model building
disc_data_id = "20260501_1246";
disc_data_folder = "/Users/leonie.thomas/Documents/fastcore_workflow_with_vanille/discretization" + filesep +  disc_data_id + filesep;
disc_data_obj = disc_data_folder + disc_data_id + "_disc_data.mat";

def_run_file = disc_data_folder + disc_data_id + "_def_run_paramters.txt";

% load needed libraries to run script
addpath(genpath("/Users/leonie.thomas/Documents/fastcore_workflow_with_vanille"))
addpath(genpath("C:\Users\leonie.thomas\rFASTCORMICS"))
addpath(genpath("C:\Users\leonie.thomas\scFASTCORMICS"))
%initCobraToolbox -> init cobratoolbox if needed
changeCobraSolver("gurobi")

%% read in all the script parameters and set working directory, directory the discretization is saved into

exp = fastcore_experiment(def_run_file); % create fastcore experiment to store models in 
load(disc_data_obj) % load data to constrain the models with

addpath(genpath(exp.script_parameters.set_working_directory));
cd(exp.script_parameters.set_working_directory);

% define a unique identifier to be used to name the different  discretizations created
date = char(datetime('now', 'Format', 'yyyyMMdd_hhss')); % to name the model and all the output  

%% start a fastcore experiment & test the consistency of input model

model = load(exp.script_parameters.model_used); %load reference model to use
model = model.(string(fieldnames(model))); % store read in model under the variable name model
model = generateRules(model); % generate model.rules field


% BUILD generic CONSISTENT model - fast consistency check (fastcc)
exp = exp.add_model(model) % 
clear model

%% add disc data 

exp = exp.add_expression_data(data,disc_data_obj)
clear data disc_data_obj

%% create medium for the models % do it before 

exp = exp.add_medium(["DMEM.tsv", "FBS.tsv"],"./data/media");

exp = exp.read_medium_files(["Concentration_M","Concentration_M"]); 
exp.medium.medium_composition.Flux_mmol_gDW_h = - conc2Rate(exp.medium.medium_composition.Concentration_mM,...
                                                       exp.script_parameters.cell_conc_cells_per_ml,...
                                                       exp.script_parameters.t_in_hours,...
                                                       exp.script_parameters.cell_dry_weight_gDW); % TODO put the values from the description file here 

 exp = exp.add_additional_rxns_boundaries(string(split(exp.script_parameters.unwanted_uptakes_export_lb, ";"))',...
                                         string(split(exp.script_parameters.unwanted_uptakes_export_ub, ";"))',...
                                         string(split(exp.script_parameters.needed_mets_medium, ";"))',...
                                         [])

exp = exp.medium_constrain("set_fluxes",0)
clear med

% group the non model components into one settings struct 

%% 
%%%%%%%%%%%%%%%%%%%%%%  set settings used for the fastcormics run 
% the unpenalized argument is currently not part of fastcormics4cobra_v2
%optional_settings.unpenalized = model_orig.rxns(ismember(vertcat(model_orig.subSystems{:}), ... 
%                                                         strsplit(exp.script_parameters.unpenalizedSystems,";")));

% forcing the medium in by setting it into fun option
biomass_rxn = exp.script_parameters.objective_function;
optional_settings.func = {'DM_atp_c_', biomass_rxn{:}, exp.medium.medium_composition.ExRxns_Recon3D{:}}; % depends if you want to foce in the medium, even though the model does not need it 
% depends on what you want how confident you are in your uptakes -> are
% they actually forced in ? not ? 
optional_settings.medium = exp.medium.medium_composition.Mets_Recon3D;
optional_settings.not_medium_constrained = exp.script_parameters.not_medium_constrained;

%%%%%%%%%%%%%%%%%%%%%%  
                        
condition_models = struct();
condition_models.orig_model.model = exp.orig_model;
condition_models.orig_model.settings.dico = exp.dico;
condition_models.consistent_medium_constrained_model.model = exp.consistent_medium_constrained_model;
condition_models.consistent_medium_constrained_model.settings.dico = exp.dico;

condition_column = exp.script_parameters.columns_to_define_model_samples_on;
% get the index of the samples in every defined group
for cond = unique(exp.data.metadata.(condition_column))'
         % transform the array, the for loop loops over the rows, so if the elements over which you want to loop over are defined in cells in one row /not column then the for loop will concat all elements instead of looping over them
        idx = contains(exp.data.metadata.(condition_column),cond);
        
        % write diary
        fname = [tempname '.txt'];
        diary(fname)
        diary on

        disp("condition for which the samples are filtered: " + cond + newline + " ----------------####################### ------------------------");
        model_cond = struct();
        % run rfastcormics on consistent global metabolic model
        tic; % mearuse the time the model takes to run
        model_cond.expression_data = table(string(exp.data.norm_counts.Properties.RowNames),...
                                            exp.data.norm_counts{:,idx},...
                                            'VariableNames',["gene_names","values"]);% add the data used for the model to the resulting model
        model_cond.discretized_data = table(string(exp.data.norm_counts.Properties.RowNames),...
                                            exp.data.discretized(:,idx),...
                                            'VariableNames',["gene_names","values"]);
        model_cond.sample_metadata = exp.data.metadata(idx,:); % add metadta of the samples used to compute the model!s
        reference_model = "consistent_medium_constrained_model";
        [model_cond.model,~, core_reaction_indices] = rFastcormics4cobra_v2(exp.(reference_model),model_cond.discretized_data.values, ...
                                                                                    model_cond.discretized_data.gene_names, exp.dico,...
                                                                                    exp.script_parameters.consensus_proportion, ...
                                                                                    exp.script_parameters.epsilon,...
                                                                                    optional_settings, biomass_rxn, 0, 0);
        model_cond.core_reactions = string(exp.(reference_model).rxns(core_reaction_indices));
        model_cond.settings = struct();
        
        %% put dico in the same order as model.genes slot 

        [isFound, idx_dico] = ismember(regexprep(string(model_cond.model.genes), "\.[0-9]+$", ""), ...
                                  regexprep(string(exp.dico.ENTREZ), "\.[0-9]+$", ""));
        dico_aligned = array2table(repmat("NA", numel(isFound), width(exp.dico)), ...
                                   'VariableNames', exp.dico.Properties.VariableNames);
        dico_aligned(isFound,:) = varfun(@string, exp.dico(idx_dico(isFound), :));
        dico_aligned.gene_id_in_model = string(model_cond.model.genes);

        %%%%% add additional settings information
        model_cond.settings.dico = dico_aligned;
        model_cond.settings.medium = exp.medium;
        model_cond.settings.obj_function = string(biomass_rxn);
        model_cond.settings.optional_settings = optional_settings;
        model_cond.settings.script_parameters = exp.script_parameters;
        model_cond.settings.reference_model = reference_model;
        model_cond.settings.running_time = toc;

        if ~sum(contains(model_cond.model.rxns,"biomass_reaction")) == 1
            error("Biomass_rxn is lost!!")
        end
        if length(model_cond.model.rxns)~=length(fastcc(model_cond.model, exp.script_parameters.epsilon, 0, 0, 'original'))
            error("Model generated by fascormics is not consistent!")
        end
        disp("number of samples for which this condition was modelled on: " + size(exp.data.discretized(:,idx),2) + newline + " ----------------####################### ------------------------")
        diary off

        model_cond.settings.diary = fileread(fname);
        model_cond.analysis = struct();
        condition_models.(strrep(cond{:},"-","_")) = model_cond;
        clear model_cond
        delete(fname);
end

exp.optional_settings = optional_settings;

clear biomass_rxn

%%

clear idx model_cond condition_column cond optional_settings fname cond
cd(exp.script_parameters.set_working_directory)
mkdir(fullfile(exp.script_parameters.save_models_to, date))

project = struct();
project.models = condition_models;

project.comparisons = struct();


project_file_name = fullfile(exp.script_parameters.save_models_to,date,[ date '_project.mat']);  % Convert datetime object to string
save(project_file_name, 'project')

exp_file_name = fullfile(exp.script_parameters.save_models_to,date,[ date '_fastcore_exp.mat']);  % Convert datetime object to string
disp(exp_file_name);
md_file_name = fullfile(exp.script_parameters.save_models_to,date,[ date '_cond_models.mat']);  % Convert datetime object to string
disp(md_file_name);
dat_file_name = fullfile(exp.script_parameters.save_models_to,date,[ date  '_workspace_cond_models.mat']);  % Convert datetime object to string
disp(dat_file_name);

save(md_file_name, 'condition_models')
clear condition_models
save(dat_file_name)
save(exp_file_name, 'exp')

[~, name, ext] = fileparts(def_run_file);
def_file_name = fullfile(exp.script_parameters.save_models_to,date, date + "_" + name + ext);
copyfile(def_run_file, def_file_name);
