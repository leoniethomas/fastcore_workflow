%% Normalization and QC of the gene expression data + Discretization 


%% set up session & define input parameter file

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

% read in the parameters needed for the analysis
%def_run_file = "\\atlas.uni.lux\fstc_sysbio\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model\data\def_run_paramters.txt";
def_run_file = "/Users/leonie.thomas/Documents/fastcore_workflow_with_vanille/data/def_run_paramters.txt";
addpath(genpath("/Users/leonie.thomas/Documents/fastcore_workflow_with_vanille"));
addpath(genpath("/Users/leonie.thomas/scFASTCORMICS"));


%% read in all the script parameters and set working directory, directory the discretization is saved into

scr_para = read_in_run_def_file(def_run_file);

% add the working path to the path & set the github repo location to be the working directory
addpath(genpath(scr_para.set_working_directory));
cd(scr_para.set_working_directory);

% define a unique identifier to be used to name the different  discretizations created
date = char(datetime('now', 'Format', 'yyyyMMdd_hhss')); % to name the model and all the output 
    
 
model = load(scr_para.model_used);
model = model.(string(fieldnames(model)));
model = generateRules(model);
model.id = scr_para.model_used;

load(scr_para.gene_dic_file)

%% load metadata, and expression data

% - as an input for fastcormics RNAseq data is used
% - normalized counts (for example vst normalized) can be used for the
%   visualization, as an input for fastcormics fpkm/tpm can be used 
% - for both the normalization is done per sample between genes, therfore
%   the visualzation with PCA/UMAP might not recover the expected cluster
%   structure 

data = expression_data(scr_para.count_data,... % the raw counts
                       scr_para.expression_data_metadata_file, ... % the metadata characterizing the samples, if available
                       dico,...
                       "sample"); % column defining the samle names 
                   
%% discretize expression data

data = data.get_normalized_data(scr_para.expression_data_expr_file);

%%

data = data.get_discretized_data(1,...
                                 string(scr_para.save_disc_data_to) + date + filesep,...
                                 '.svg');
                             

%% mapping gene expression to rxns 

data = data.map_expression_2_rxns(model,dico)

%% load model and define metabolic genes for visualization 


data = data.get_metabolic_genes(load(scr_para.model_used),...
                                "SYMBOL");

%% 

% creating a directory for the following discretzation to be stored into - copying the defined parameters 
mkdir((scr_para.save_disc_data_to), date)
copyfile(def_run_file, ...
         string(scr_para.save_disc_data_to) + date+ filesep + date + "_def_run_paramters.txt")
     

clear def_run_file input_paramters 


%% save data for model construction
dat_file_name = fullfile(scr_para.save_disc_data_to,   date , [date  '_disc_data.mat']);  % Convert datetime object to string
disp(dat_file_name);

save(dat_file_name, 'data')

%% export data object ? so it can be read in by R ? - see the R function on my spatial models from  brussels ? 
outdir = scr_para.save_disc_data_to + filesep + date + filesep +  "exported_csv";
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

save_struct_to_folder(struct(data),outdir)


% %% QC - data exploration
% 
% 
% % run pca on data 
% data = data.perform_pca_kmeans("FPKM",2);
% data = data.perform_pca_kmeans("mapping_exp_2_rxns",2);
% data = data.perform_pca_kmeans("TPM",2);
% data = data.perform_pca_kmeans("discretized",2);
% 
% %% perform umap and visualize the umap with visualize_dimreduction function
% 
% addpath(genpath(string(scr_para.set_working_directory) + filesep + "scr"+filesep + "fastcormics_module" +filesep + "UMAP"));
% rehash;
% %data = data.perform_umap("FPKM")
% 
% %data.visualize_dimreduction("kmeans_k2_FPKM_features_37232","umap","FPKM")
% 
% 
% %%
% % visualize pca and clustering
% data.visualize_dimreduction("kmeans_k2_FPKM_features_37232","pca","FPKM")
% data.visualize_dimreduction("Treatment","pca","FPKM")
% %data.visualize_dimreduction("kmeans_k2_mapping_exp_2_rxns_features_10600","pca","mapping_exp_2_rxns")
% data.visualize_dimreduction("Treatment","pca","mapping_exp_2_rxns")
% data.visualize_dimreduction("kmeans_k2_TPM_features_37232","pca","TPM")
% data.visualize_dimreduction("Treatment","pca","TPM")
% data.visualize_dimreduction("kmeans_k2_discretized_features_37232","pca","discretized")
% data.visualize_dimreduction("Treatment","pca","discretized")
% 
% 
% 
% % get QC plots                               
% data.get_QC_plots("raw_counts","SampleID",string(scr_para.QC_figures_path) + date + filesep)
% 
% data.get_QC_plots("FPKM","SampleID",string(scr_para.QC_figures_path) + date + filesep)
% 
