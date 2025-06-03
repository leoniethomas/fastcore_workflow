%% Normalization and QC of the gene expression data used as an input to constrain the metabolic model and therefore create a context specific model

% - this scripts works with an file defining the input parameters - downloading the github repository the file can be found in the data folder 


%% read in all the script parameters and set working directory, directory the discretization is saved into

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

% read in the parameters needed for the analysis
def_run_file = "\\atlas.uni.lux\fstc_sysbio\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model\data\def_run_paramters.txt";
input_paramters = readtable(def_run_file, 'Delimiter','\t');
scr_para = cell2struct(input_paramters{:,"value"}, input_paramters{:,"slot_name"});

% add the working path to the path & set the github repo location to be the
% working directory
addpath(genpath(scr_para.set_working_directory));
cd(scr_para.set_working_directory);

% define a unique identifier to be used to name the different
% discretizations created
date = char(datetime('now', 'Format', 'yyyyMMdd_hhss')); % to name the model and all the output 

% creating a directory for the following discretzation to be stored into -
% copying the defined parameters 
mkdir((scr_para.save_disc_data_to), date)
mkdir((scr_para.QC_figures_path), date)
copyfile(def_run_file, ...
         string(scr_para.save_disc_data_to) + date+ filesep + date + "_def_run_paramters.txt")
     

clear def_run_file input_paramters 



%% load metadata, and expression data

% - as an input for fastcormics RNAseq data is used
% - normalized counts (for example vst normalized) can be used for the
%   visualization, as an input for fastcormics fpkm/tpm can be used 
% - for both the normalization is done per sample between genes, therfore
%   the visualzation with PCA/UMAP might not recover the expected cluster
%   structure 

data = expression_data(scr_para.count_data,... % the raw counts
                       scr_para.expression_data_metadata_file, ... % the metadata characterizing the samples, if available
                       "SampleID"); % column defining the samle names 
                   
%% discretize expression data

data = data.get_normalized_data(scr_para.expression_data_expr_file,"FPKM");
data = data.get_normalized_data(scr_para.TPM_file,"TPM");

%%

data = data.get_discretized_data(1,...
                                 string(scr_para.QC_figures_path) + date + filesep, ...
                                 "TPM");
                             
%% 

% - visualize with only metabolic genes -> 

%% load model and define metabolic genes for visualization 


data = data.get_metabolic_genes(load(scr_para.model_used),...
                                load(scr_para.gene_dic_file));

%% QC - data exploration


[coeff,score,latent,...
 tsquared,explained,cluster] = data.QC_pca_kmeans("TPM", scr_para.columns_to_define_model_samples_on, ...
                                                  [1 2],6, 1:length(data.feature_names_norm), ...
                                                  [scr_para.QC_figures_path date '/' date  'TPM_PCA.png'],...
                                                  ["MDA-MB231-" ""]);
[coeff,score,latent,...
 tsquared,explained,cluster] = data.QC_pca_kmeans("discretized", scr_para.columns_to_define_model_samples_on, ...
                                                  [1 2],2,1:length(data.feature_names_norm), ...
                                                  [scr_para.QC_figures_path date '/' date  '_discretized_on_TPM_metgenes_PCA.png'],...
                                                  ["MDA-MB231-" ""]);   
                                              
[coeff,score,latent,...
 tsquared,explained,cluster] = data.QC_pca_kmeans("FPKM", scr_para.columns_to_define_model_samples_on, ...
                                                  [1 2],6, 1:length(data.feature_names_norm), ...
                                                  [scr_para.QC_figures_path date '/' date  '_fpkm_PCA.png'],...
                                                  ["MDA-MB231-" ""]);   
                                              
[coeff,score,latent,...
 tsquared,explained,cluster] = data.QC_pca_kmeans("raw_counts", scr_para.columns_to_define_model_samples_on, ...
                                                  [1 2],6, 1:length(data.feature_names_raw), ...
                                                  [scr_para.QC_figures_path date '/' date  '_raw_counts_PCA.png'],...
                                                  ["MDA-MB231-" ""]);                                          


data.get_QC_plots("raw_counts","SampleID",string(scr_para.QC_figures_path) + date + filesep)

data.get_QC_plots("FPKM","SampleID",string(scr_para.QC_figures_path) + date + filesep)

%% save data for model construction
dat_file_name = [scr_para.save_disc_data_to   date '/' date  '_disc_data.mat'];  % Convert datetime object to string
disp(dat_file_name);

save(dat_file_name, 'data')


