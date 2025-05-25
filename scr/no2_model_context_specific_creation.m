%% analysis of TNBC metabolism


%% Setup and parameters


clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

% read in the parameters needed for the analysis
vis_only_met_genes = 1;
def_run_file = "\\atlas.uni.lux\fstc_sysbio\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model\data\def_run_paramters.txt";
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

%% load metadata, and expression data

data = expression_data(scr_para.count_data,... % the raw counts
                       scr_para.expression_data_metadata_file, ... % the metadata characterizing the samples, if available
                       "SampleID"); % column defining the samle names 
                   
%% discretize expression data
data = data.get_normalized_data(scr_para.expression_data_expr_file); 
%%
data = data.get_discretized_data(1,...
                                 string(scr_para.QC_figures_path) + date + filesep);

%% load model and define metabolic genes for visualization 

load(scr_para.model_used) 
load(scr_para.gene_dic_file)
metabolic_genes_entrez = string(regexprep(model.genes,".1",""));
metabolic_genes_mapped_dico = string(dico.SYMBOL(find(matches(dico.ENTREZ,metabolic_genes_entrez))));
metabolic_gene_idx = find(matches(metabolic_genes_mapped_dico, data.feature_names_norm));


%% QC - data exploration


[coeff,score,latent,...
 tsquared,explained,cluster] = data.QC_pca_kmeans("normalized_counts", scr_para.columns_to_define_model_samples_on, ...
                                                  [1 2],6, metabolic_gene_idx, ...
                                                  [scr_para.QC_figures_path date '/' date  'normalized_counts_PCA_.png'],...
                                                  ["MDA-MB231-" ""]);
                                              
[coeff,score,latent,...
 tsquared,explained,cluster] = data.QC_pca_kmeans("raw_counts", scr_para.columns_to_define_model_samples_on, ...
                                                  [1 2],6, 1:length(data.feature_names_raw), ...
                                                  [scr_para.QC_figures_path date '/' date  '_raw_counts_PCA_.png'],...
                                                  ["MDA-MB231-" ""]);                                          


data.get_QC_plots("raw_counts","SampleID",string(scr_para.QC_figures_path) + date + filesep)

data.get_QC_plots("normalized_counts","SampleID",string(scr_para.QC_figures_path) + date + filesep)



%% save all the figures to a separate folder 
%copyfile('Figures/QC', ['Figures/QC/' date])
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
%% DISCRETIZE expresssion data 

expression_data.discretized = discretize_FPKM(expression_data.data, expression_data.sample_names,1);
% movefile('Figures/Discretization', ['Figures/QC/Discretization/' date])

%discretized = discretize_FPKM(expression_data.data, expression_data.colnames,1); %with figures, will save figures in Figures folder


clear num_disc perc_disc 

%% check the discretized data for the grouping

[coeff,score,latent,tsquared,explained] = pca(expression_data.discretized(vis_genes_idx,:)');

count_data.metadata.cluster =num2cell(num2str(kmeans(expression_data.discretized(vis_genes_idx,:)',6)));


[count_data.metadata.sil,~] = silhouette(expression_data.discretized(vis_genes_idx,:)',...
                                         count_data.metadata.cluster,'Euclidean');
                                     
disp("###----- silhoutte value: " + num2str(mean(count_data.metadata.sil)) + " -----###")

disp("###----- number of samples sil <0 : " + num2str(sum(count_data.metadata.sil<0)) + " -----###") 

if mean(count_data.metadata.sil)>0.1 
    %count_data.metadata.cluster_dis = regexprep(count_data.metadata.Groups, "MDA-MB231-", "") + "_" + count_data.metadata.cluster_dis;

    figure
    hold on
    for x = unique(count_data.metadata.cluster)'
        disp(x);
        idx = contains(count_data.metadata.cluster,x{:});
        scatter(score(idx,1),score(idx,2))
        %text(score(idx,1),score(idx,2), count_data.metadata.SampleID(idx), ...
        %             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end   
    title('PCA on discretized data - fpkm')

    legend(unique(count_data.metadata.cluster)' ,'location','best')
else
    disp("###----- silhoutte value is below 0.1 - therefore the clusters are overlapping, does not really make sense to visualize them! ----###")
end


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
for cond = unique(expression_data.metadata.(condition_column))'
         % transform the array, the for loop loops over the rows, so if the elements over which you want to loop over are defined in cells in one row /not column then the for loop will concat all elements instead of looping over them
        idx = contains(expression_data.metadata.(condition_column),cond);

        disp("condition for which the samples are filtered: " + cond + newline + " ----------------####################### ------------------------");
        
        % run rfastcormics on consistent global metabolic model
        tic; % mearuse the time the model takes to run
        [model_cond,AA] = fastcormics_RNAseq(model_orig,expression_data.discretized(:,idx), ...
                                             expression_data.feature_name, dico, biomass_rxn, str2double(scr_para.already_mapped_tag),...
                                                str2double(scr_para.consensus_proportion), str2double(scr_para.epsilon), optional_settings);
        model_cond.running_time = toc;
        model_cond.used_data = expression_data.discretized(:,idx); % add the data used for the model to the resulting model
        model_cond.sample_metadata = expression_data.metadata(idx,:); % add metadta of the samples used to compute the model!
        model_cond.AA = AA;
        condition_models.(strrep(cond{:},"-","_")) = model_cond;
        disp("number of samples for which this condition was modelled on: " + size(expression_data.discretized(:,idx),2) + newline + " ----------------####################### ------------------------")
end


%%

clear A AA idx xi x TXT tsquared model_cond condition_column cond

md_file_name = [scr_para.save_models_to  date '/' date '_cond_models.mat'];  % Convert datetime object to string
disp(md_file_name);
dat_file_name = [scr_para.save_models_to   date '/' date  '_workspace_cond_models.mat'];  % Convert datetime object to string
disp(dat_file_name);

save(md_file_name, 'condition_models')
save(dat_file_name)


