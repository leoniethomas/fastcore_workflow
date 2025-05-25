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
mkdir((scr_para.save_disc_data_to), date)
mkdir((scr_para.QC_figures_path), date)

copyfile(def_run_file, ...
         string(scr_para.save_disc_data_to) + date+ filesep + date + "_def_run_paramters.txt")
     

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

%% save data for model construction
dat_file_name = [scr_para.save_disc_data_to   date '/' date  '_disc_data.mat'];  % Convert datetime object to string
disp(dat_file_name);

save(dat_file_name, 'data')

