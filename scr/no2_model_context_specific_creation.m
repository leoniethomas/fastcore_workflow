%% analysis of NO metabolism in bulk RNA seq data

% - publically available dataset (see:
% https://www.nature.com/articles/s41588-021-00911-1, GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz)
  

% This script is based on the melanoma running example on the sysbio github page (see https://github.com/sysbiolux/rFASTCORMICS/tree/master/rFASTCORMICS%20for%20RNA-seq%20data/RunningExample_Melanoma)
% as well as the DataConversionAndQC script -> https://github.com/sysbiolux/ISB705MetabolicNetworkModeling/blob/main/DataProcessing/MATLAB/DataConversionAndQC.m
%       - this DataConversion script also entails the conversion from counts ->  pseudocounts -> FPKM or TPM
%       - here the count data was processed in R with script no1, so the expression data is directly read in as FPKMs

%% ToDo: 
% - [x] create a good README file in every folder!!!
% - [ ] installation instruction!!! information about the matlab and cplex version used
% - [ ] write into the README 
% - [ ] rename the files


% normalization data -> “LogNormalize”: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p

%% Setup and parameters


clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

disp("###----- Setting all variables needed in the script -----###") 

set_working_directory = 'C:\Users\leonie.thomas\20241104_TNBC_wu_2021';
cd(set_working_directory);
save_models_to = "./models/";
QC_figures_path = "./Figures/QC/";

% model and medium (see:
% https://github.com/sysbiolux/rFASTCORMICS/tree/master/rFASTCORMICS%20for%20RNA-seq%20data/RunningExample_Melanoma)
model_used="./data/Recon3DModel.mat"; % define generic input model used
% 
model_description='Recon';

medium_used.file="./data/RPMI_Formulation.xlsx";
medium_used.concentration_column='mM';
medium_used.naming_colum='Components';


% question: why are exactly those unwanted ? 
unwanted_uptakes_export.ub = ["EX_oh1[e]", ... 
                              "sink_band[c]"];
unwanted_uptakes_export.lb = ["EX_h2o2[e]", ...
                              "EX_o2s[e]", ...
                              "EX_oh1[e]", ...
                              "EX_ppi[e]", ...
                              "sink_fe3[c]", ...
                              "sink_band[c]"];

% a few of the reactions are closed because the reactions are not realistic
% for example geting the oxide from competitiors like o2s or ohl 
% therefore those are closed 
% intracellular iron production is also unrealistic which is why the cells
% should take it from the medium, therefore the sink reaction fe is closed 
                       
                          
% expression data parameters 
expression_data.expr_file = './data/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts_fpkm.csv';%
%expression_data.expr_file = './data/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts_normalized.csv';%
%expression_data.expr_file = './data/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts_normalized_DeSeq.csv';%
%expression_data.expr_file = './data/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt';
expression_data.metadata_file = './data/GSE176078_Wu_etal_2021_bulkRNAseq_metadata.txt';
count_data = './data/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt';

% dictionary defining symbol to entrezid gene expression 

gene_dic_file = './data/dicorFASTCORMICS.mat';

% parameters context specific model generation

epsilon = 1e-4;
consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included.
% Only relevant if you want to create one model from different samples
already_mapped_tag = 0;

unpenalizedSystems = {'Transport, endoplasmic reticular'; % Question: what does unpeanalized mean ??
                      'Transport, extracellular';
                      'Transport, golgi apparatus';
                      'Transport, mitochondrial';
                      'Transport, peroxisomal'; % question: how is this selected ? 
                      'Transport, lysosomal';
                      'Transport, nuclear'};

                  
% high density lipo protein - this is one of the lipids that is missing in
% the model, which the cell needs for the lipid metabolism to work properly
% therefore it is added from external
not_medium_constrained = {'EX_hdl_hs[e]'};


columns_to_define_model_samples_on = ["condition"]% ,"CaseID"]; 

date = char(datetime('now', 'Format', 'MMddyyyy')) % to name the model and all the output 

%% load metadata
disp("###----- Load metadata -----###") 
sample_metadata = readtable(expression_data.metadata_file, 'Delimiter','\t');

sample_metadata.CaseID = strcat("CID",replace(sample_metadata.CaseID,"-",""));
sample_metadata.new_sample_id = strcat(sample_metadata.CaseID,"_",sample_metadata.SubtypeByIHC);
sample_metadata.condition = sample_metadata.SubtypeByIHC;

for condition_column= columns_to_define_model_samples_on
    sample_metadata.(condition_column) =  regexprep(sample_metadata.(condition_column),'[^a-zA-Z0-9]', '');
end

%% load expression data and perform QC

disp("###----- Load expression data -----###") 
% load FPKM data

expression_data.data = readcell(expression_data.expr_file); % read in as string to also get colnames

% get colnames from the read in dataframe
colnames = ["ID", ...
            expression_data.data(1,...
                                 cellfun(@(x) ~isempty(regexp(x, '[A-Za-z]')), string(expression_data.data(1,:)), ...
                                 'UniformOutput', ...
                                 true)...
                                )...
           ]; % get colnames from first row


disp("###----- align expression data with sample metadata -----###") 
% bring the samples in the same order in the expression as well as the metadata
[~,A,B] = intersect(colnames,sample_metadata.CaseID);
colnames = colnames([1,A']); % add the first column because gene id 
expression_data.sample_metadata = sample_metadata(B,:);
expression_data.data = expression_data.data(:,[1,A']); % add first column, which is the gene id
expression_data.data(1,:) = []; % remove first column, is later stored in colnames 
expression_data.data = cell2table(expression_data.data, 'VariableNames',colnames); % transform the expression data to table, only numbers, no col or rownames

% save colnames and rownames to the expression data object
expression_data.rownames=expression_data.data.ID;
rownames=strtok(expression_data.rownames,'.');
expression_data.data = table2array(expression_data.data(:,2:length(colnames)));
expression_data.colnames = colnames(2:length(colnames));

%clear A B sample_metadata colnames rownames

%% load count data

counts = readtable(count_data);

counts = table2array(counts(:,2:end)); %log2 scaling with pseudocount 1
counts = counts(:,A-1');
disp('Number of reads:')
sum(counts)
counts = log2(counts+1); %log2 scalling with pseudocount 1
mapcaplot(counts',expression_data.colnames);


%% perform QC

disp("###----- perform QC -----###") 
disp("###----- PCA")

% with pca function and manual plotting
[~,score,~,~,explained,~] = pca(expression_data.data');
pc1 = 1;
pc2 = 2;
catLabels = unique(expression_data.sample_metadata.(condition_column)');
colors = {"g", "b", "k", "r", "b"};

figure
hold on
for x = catLabels
    x = string(x);
    idx = find(expression_data.sample_metadata.condition'==x);
    dot_col = colors{find(x == catLabels)};
    plot(score(idx,pc1),score(idx,pc2), dot_col + '.', 'MarkerSize',20 )
end   

title('PCA')
xlabel([num2str(pc1), ' component: ', num2str(explained(pc1))])
ylabel([num2str(pc2), ' component: ', num2str(explained(pc2))])
legend(catLabels ,'location','best')

saveas(gcf, QC_figures_path + date + '_fpkm_PCA_.png');

% with pca function and manual plotting
[~,score,~,~,explained,~] = pca(counts');
pc1 = 1;
pc2 = 2;
catLabels = unique(expression_data.sample_metadata.(condition_column)');
colors = {"g", "b", "k", "r", "b"};

figure
hold on
for x = catLabels
    x = string(x);
    idx = find(expression_data.sample_metadata.condition'==x);
    dot_col = colors{find(x == catLabels)};
    plot(score(idx,pc1),score(idx,pc2), dot_col + '.', 'MarkerSize',20 )
end   

title('PCA')
xlabel([num2str(pc1), ' component: ', num2str(explained(pc1))])
ylabel([num2str(pc2), ' component: ', num2str(explained(pc2))])
legend(catLabels ,'location','best')

saveas(gcf, QC_figures_path + date + '_PCA_.png');

disp("###----- barplot expression per sample")
bar(sum(expression_data.data))
title('Number of reads per sample: ')
xlabel("samples")
ylabel("# of reads")
saveas(gcf, QC_figures_path + date + '_expression_read_counts_barplot.png');

bar(sum(expression_data.data == 0,1))
title('Number of 0 per sample')
xlabel("samples")
ylabel("# of reads")
saveas(gcf, QC_figures_path + date + '_expression_zero_counts_barplot.png');

disp("###----- boxplot expression per sample")

figure
boxplot(counts)
saveas(gcf, QC_figures_path + date + '_expression_unnormalized_boxplots.png');

figure
boxplot(expression_data.data)
saveas(gcf, QC_figures_path + date + '_expression_normalized_boxplots.png');

disp("###----- boxplot expression per sample - zeros")

% without zeros
data2=expression_data.data;
data2(data2==0)=NaN;
figure
boxplot(data2)

saveas(gcf, QC_figures_path + date + '_expression_normalized_boxplots_nozeros.png');

clear data2

% without zeros
data2=counts;
data2(data2==0)=NaN;
figure
boxplot(data2)

saveas(gcf, QC_figures_path + date + '_expression_unnormalized_boxplots_nozeros.png');

clear data2


disp("###----- density plots of normalized data per sample")

figure
hold on
fpkm = counts;
fpkm(fpkm==0)=NaN; %remove zeros for densityplot
fpkm=log2(fpkm); %log2 scaling
for i=1:size(fpkm,2)
    [probability_estimate,xi] = ksdensity(fpkm(:,i));
    plot(xi,probability_estimate,':k','LineWidth',1);
end
hold off
saveas(gcf, QC_figures_path + date + '_density_dist_unnormalized_data.png');

figure
hold on
fpkm = expression_data.data;
fpkm(fpkm==0)=NaN; %remove zeros for densityplot
fpkm=log2(fpkm); %log2 scaling
for i=1:size(fpkm,2)
    [probability_estimate,xi] = ksdensity(fpkm(:,i));
    plot(xi,probability_estimate,':k','LineWidth',1);
end
hold off
saveas(gcf, QC_figures_path + date + '_density_dist_normalized_data.png');


%% Load medium

[NUM,TXT,RAW]=xlsread(medium_used.file);
medium_RPMI=TXT(2:end,2);
medium_RPMI_EX=cellfun(@(x)['EX_' x],medium_RPMI,'uni',false);
T=table(medium_RPMI_EX,medium_RPMI,NUM(:,3));


%% LOAD MODEL: Recon and remove unwanted export and take up reactions



load(model_used) % where does this file come from ? 
model.description=model_description;


t = readtable("./data/fluxes.tsv", "FileType","text",'Delimiter', '\t');
t = t(find(t.FluxValue < 1000),:) % does not really make sense if it is higher
t.FluxValue = - t.FluxValue;

% find the intersection of reactions in the flux file, model & the medium
fluxes_to_constrain_in_model = intersect(model.rxns(findExcRxns(model)),...
                                                   t.Reaction)
%fluxes_to_constrain_in_model = intersect(intersect(model.rxns(findExcRxns(model)),...
%                                                   t.Reaction),...
%                                         T.medium_RPMI_EX)
[~,idx, idx_fluxes_in_model] = intersect(fluxes_to_constrain_in_model,model.rxns);
[~,idx, idx_fluxes_in_flux_file] = intersect(fluxes_to_constrain_in_model,t.Reaction);                                   

model.lb(idx_fluxes_in_model) = t{idx_fluxes_in_flux_file,"FluxValue"};


model.ub(find(ismember(model.rxns,unwanted_uptakes_export.ub)))=0; 
model.lb(find(ismember(model.rxns,unwanted_uptakes_export.lb)))=0; 



%% BUILD generic CONSISTENT model - fast consistency check (fastcc)

A = fastcc_4_rfastcormics(model, 1e-4, 1);
% or with COBRA fastcc:
% [A, ~, ~] = fastcc(model, 1e-4, 1, 0, 'original');

% remove non consistent reactions from model
model=removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));
% check if the biomass reactions are still there
model.rxns(find(contains(model.rxns,'biomass')))

model_orig=model;
% check if the created model is now really consistent
A = fastcc_4_rfastcormics(model, 1e-4, 1);

clear model


%% DISCRETIZE expresssion data 

expression_data.discretized = discretize_FPKM(expression_data.data, expression_data.colnames,1); % in cases where the left distribution is higher than the right one
%discretized = discretize_FPKM(expression_data.data, expression_data.colnames,1); %with figures, will save figures in Figures folder

figure
hist(expression_data.discretized)
tabulate(expression_data.discretized(:,1))

%% BUILT CONTEXT SPECIFIC MODELS -> reconstruction using rFASTCORMICS

load(gene_dic_file)
subSys=vertcat(model_orig.subSystems{:});
 
optional_settings.unpenalized = model_orig.rxns(ismember(subSys,unpenalizedSystems));
optional_settings.func = {'DM_atp_c_', 'biomass_reaction', T.medium_RPMI_EX{:}}; %biomass_maintenance %-> c
optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = T.medium_RPMI; %(add media instead)
biomass_rxn = {'biomass_reaction'} 

                        
condition_models = struct();

for condition_column= columns_to_define_model_samples_on % looping over all columns which were entail categories deviding the samples into groups, used to compute consensus models of the defined group 
    % here condition so disease condition and the indiviudal sample
    % indentifier were choosen
    idx_cond = arrayfun(@(x) find(strcmp(expression_data.sample_metadata.(condition_column), x)), ...
                        unique(expression_data.sample_metadata.(condition_column)),'UniformOutput',false); 
    % get the index of the samples in every defined group
    for idx = idx_cond' % transform the array, the for loop loops over the rows, so if the elements over which you want to loop over are defined in cells in one row /not column then the for loop will concat all elements instead of looping over them
        idx = cell2mat(idx)'; % transform so it can be used to index 
        cond = cell2mat(unique(expression_data.sample_metadata(idx,:).(condition_column))); % what is the condition of the indexed samples? 
        disp("condition for which the samples are filtered: " + cond + newline + " ----------------####################### ------------------------");
        
        % run rfastcormics on consistent global metabolic model
        tic; % mearuse the time the model takes to run
        [model_cond,AA] = fastcormics_RNAseq(model_orig,expression_data.discretized(:,idx), ...
                                             expression_data.rownames, dico, biomass_rxn, already_mapped_tag,...
                                             consensus_proportion, epsilon, optional_settings);
        model_cond.running_time = toc;
        model_cond.used_data = expression_data.discretized(:,idx); % add the data used for the model to the resulting model
        model_cond.sample_metadata = expression_data.sample_metadata(idx,:); % add metadta of the samples used to compute the model!
        model_cond.AA = AA;
        condition_models.(cond) = model_cond;
        disp("number of samples for which this condition was modelled on: " + size(expression_data.discretized(:,idx),2) + newline + " ----------------####################### ------------------------")
    end
end




md_file_name = save_models_to + date + "_cond_models.mat";  % Convert datetime object to string
disp(md_file_name);
dat_file_name = save_models_to + date + "_input_data_obj_models.mat";  % Convert datetime object to string
disp(dat_file_name);
save(md_file_name, 'condition_models')
save(dat_file_name, 'expression_data')
orig_file_name = save_models_to + date + "_orig_model.mat";  % Convert datetime object to string
disp(orig_file_name);
save(orig_file_name, 'model_orig')

