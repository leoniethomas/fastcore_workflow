%% define VARIABLES 

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names



%% define script parameters
model_id = "20250227_1224";
project_path = "\\atlas.uni.lux\FSTC_SYSBIO\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model";
path_to_model_to_analyse = project_path + filesep + "context_specific_models" + filesep +model_id;
path_to_sampling_data = "." +filesep + "analysis" + filesep +model_id + filesep + "sampling" + filesep ;
cd (project_path)
addpath(genpath(project_path))
sampling_files = string(ls(path_to_sampling_data));
sampling_files = sampling_files(contains(sampling_files, "samplingR"));
sampling_files = regexprep(sampling_files, " ", "")


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
                   
                   


%% load sampling data 

sample_data = cell2struct(arrayfun(@(x) load(fullfile(x)), sampling_files','UniformOutput', false),...
                          regexprep(sampling_files, ".mat", ""),...
                          2)
models =  structfun(@(y) y.x.modelSampling,sample_data,'UniformOutput', false);
samples =  structfun(@(y) y.x.samples,sample_data,'UniformOutput', false);



samples_ordered = arrayfun(@(x) get_sampling_orig_order(models.(x),samples.(x),length(model_orig.rxns)), ...
                                  string(fieldnames(models)),...
                                  'UniformOutput',false);
                  
                              
samples_joined_ordered = cell2mat(samples_ordered');
%this helps the PCA, cause less dimensions, but in the downstream we lose
%th
% samples_joined_ordered( ~any(samples_joined_ordered,2), : ) = [];

               
%% perform PCA
visualize_sampling(string(fieldnames(models))', ...
                   regexprep(regexprep(string(fieldnames(condition_models)),"_"," "),"MDA MB231 ",""), ...
                   samples_joined_ordered',...
                   numel(string(fieldnames(models))'),...
                   2,3,...
                   1)
               
%% PCA subsystem specific 

subsystems = unique(arrayfun(@(x) x{:},model_orig.subSystems));
sil = [];
homo = []; 

for subsys_id = 1:numel(subsystems)
    % subsys = subsystems{16};
    subsys = subsystems{subsys_id}
    subsys_rxn_idx = find(string(arrayfun(@(x) x{:},model_orig.subSystems)) == subsys);
    samples_subsys_rxns =  samples_joined_ordered(subsys_rxn_idx,:);
    if length(subsys_rxn_idx) > 3
        [s,h] = visualize_sampling(string(fieldnames(models))', ...
                           regexprep(regexprep(string(fieldnames(condition_models)),"_"," "),"MDA MB231 ",""), ...
                           samples_subsys_rxns',...
                           numel(string(fieldnames(models))'),...
                           2,3,...
                           0);
    else 
        s = NaN;
        h = NaN;
    end
    sil = [sil, s];
    homo = [homo, h];              
end

%% 
subsys_id = 64;
subsys = subsystems{subsys_id};
subsys_rxn_idx = find(string(arrayfun(@(x) x{:},model_orig.subSystems)) == subsys);
samples_subsys_rxns =  samples_joined_ordered(subsys_rxn_idx,:);
[s,h] = visualize_sampling(string(fieldnames(models))', ...
                           regexprep(regexprep(string(fieldnames(condition_models)),"_"," "),"MDA MB231 ",""), ...
                           samples_subsys_rxns',...
                           numel(string(fieldnames(models))'),...
                           1,2,...
                           1);
                       
%% 
subsys_id = 24;
subsys = subsystems{subsys_id};
subsys_rxn_idx = find(string(arrayfun(@(x) x{:},model_orig.subSystems)) == subsys);
samples_subsys_rxns =  samples_joined_ordered(subsys_rxn_idx,:);
[s,h] = visualize_sampling(string(fieldnames(models))', ...
                           regexprep(regexprep(string(fieldnames(condition_models)),"_"," "),"MDA MB231 ",""), ...
                           samples_subsys_rxns',...
                           numel(string(fieldnames(models))'),...
                           1,2,...
                           1);                       
               
%% visualize the distributions for rxns

rxn_id = 2882;
rxn_samp_fluxes = cell2mat(arrayfun(@(x) x{1}(rxn_id,:),samples_ordered','UniformOutput',false)');

figure
for i=1:size(rxn_samp_fluxes,1)
    [probability_estimate,xi] = ksdensity(rxn_samp_fluxes(i,:));
    plot(xi,probability_estimate*100,'LineWidth',1); % multiplied with 100 to have %
    %trapz(probability_estimate, xi) % should approx to 1 % integral should sum up to 1
    %[y,x] = hist(rxn_samp_fluxes(i,:),25)
    %plot(x,y,'LineWidth',1);
    hold on
end
legend(regexprep(regexprep(string(fieldnames(condition_models)),"_"," "),"MDA MB231 ",""))
xlabel("rxn flux value")
ylabel("probability of obtaining x [%]")
title("Probability distributions between different models given the performed sampling - rxn: " +  model_orig.rxnNames{rxn_id} + " (idx: " + num2str(rxn_id) + " )" )
hold off



