
%% analysis of the context specific models 
% what are the main questions for now ? 

% growth rate ? is it reasonable ? 


%% define VARIABLES 

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

addpath(genpath("C:\Users\leonie.thomas\plot2svg"))

%% define script parameters
model_id = "20250525_0950";
project_path = "\\atlas.uni.lux\FSTC_SYSBIO\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model";
path_to_model_to_analyse = project_path + "\context_specific_models\" + model_id;
cd (project_path)
addpath(genpath(project_path))

%% load the created models with their whole workspace

load(path_to_model_to_analyse + "\" +   model_id + "_workspace_cond_models.mat") % load the condition specific models created with rFASTCORMICS

input_paramters = dir(path_to_model_to_analyse + "\" + "*def_run_paramters.txt");
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
                   
condition_models = rmfield(condition_models,'MDA_MB231_HERVK_C_NO')
condition_models = rmfield(condition_models,'MDA_MB231_HERVK_C_VC')
condition_models = rmfield(condition_models,'MDA_MB231_HERVK_D_NO')
condition_models = rmfield(condition_models,'MDA_MB231_HERVK_D_VC')

                   
                   
model_names = regexprep(fieldnames(condition_models),"_", " ");
%%

%writeCbModel(condition_models.MDA_MB231_Cont_NO,'format', 'json','fileName','model_Cont_NO.json')
%writeCbModel(condition_models.MDA_MB231_Cont_VC,'format', 'json','fileName','model_Cont_VC.json')

%% script initialization

changeCobraSolver("ibm_cplex");
mkdir(scr_para.results_path);
results = struct();

% load functions
fun = functions_no3;

%% 
ess_file_name = "\\atlas.uni.lux\FSTC_SYSBIO\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model\analysis\20250525_0950\ess_genes.csv"
ess_genes = readtable(ess_file_name,'ReadRowNames',true);

ess_database_file = '\\atlas.uni.lux\fstc_sysbio\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model\data\essentiall_genes\gene_effects_ACH-000768.csv';
ess_data = readtable(ess_database_file);

%% filter database for metabolic genes

met_genes = regexprep(model_orig.genes,".1$","");
[~,ia,ib] = intersect(string(met_genes), dico.ENTREZ); 
met_genes(ia) = dico.SYMBOL(ib); %extract the symbols

ess_data = ess_data(find(matches(ess_data.gene,met_genes)),:);

%% find genes which are used in the model
ctr_model = removeUnusedGenesFastbox(condition_models.MDA_MB231_Cont_VC,1);

met_genes_in_model = regexprep(ctr_model.genes,".1$","");
[~,ia,ib] = intersect(string(met_genes_in_model), dico.ENTREZ); 
met_genes_in_model(ia) = dico.SYMBOL(ib); %extract the symbols

ess_data = ess_data(find(matches(ess_data.gene,met_genes_in_model)),:);

%%

[~,ia,ib] = intersect(string(ess_genes.Properties.RowNames), dico.ENTREZ); 
ess_genes.Properties.RowNames(ia) = dico.SYMBOL(ib); %extract the symbols

%% 

ess_genes = ess_genes(find(sum(ess_genes{:,:},2)>0),:);
ess_genes = string(ess_genes.Properties.RowNames);

%%

idx_ess_genes = find(matches(ess_data.gene,ess_genes));
res_ess = ess_data(idx_ess_genes,:);

res_ess.pos_crispr_screen_met_genes = idx_ess_genes;


%% 

top20 = res_ess(res_ess.pos_crispr_screen_met_genes <= 20,:);
top20.entrez = string(dico.ENTREZ(find(matches(dico.SYMBOL,top20_symbol.gene))));

%%

ctr_model.subSystems(find(matches(ctr_model.genes,top20.entrez)))




























