%% analysis of TNBC metabolism


%% Setup and parameters


clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names

% read in the parameters needed for the analysis

def_run_file = "C:\Users\leonie.thomas\20250225_glynn_bulk_metabolic_model\data\def_run_paramters.txt";
input_paramters = readtable(def_run_file, 'Delimiter','\t');

% Define your array of field names
scr_para = cell2struct(input_paramters{:,"value"}, input_paramters{:,"slot_name"});


date = char(datetime('now', 'Format', 'MMddyyyy')); % to name the model and all the output 


%% load metadata
disp("###----- Load metadata -----###") 

cd(scr_para.set_working_directory);
sample_metadata = readtable(scr_para.expression_data_metadata_file, 'Delimiter','\t');

%% load normalized counts
disp("###----- Load expression data -----###") 

expression_data.data = readcell(scr_para.expression_data_expr_file); 
expression_data.sample_names = expression_data.data(1,:);
expression_data.feature_name = expression_data.data(2:end,1);
expression_data.data = cell2mat(expression_data.data(2:end,2:end));

% bring the metadata in the same order and then add it to object 
[~, idx, idy] = intersect(expression_data.sample_names(arrayfun(@(x) isstr(x{:}),expression_data.sample_names)),...
                          sample_metadata.SampleID);
expression_data.sample_names = expression_data.sample_names(idx);
expression_data.data = expression_data.data(:,idx);
expression_data.metadata = sample_metadata(idy,:);

clear idx idy 

%% load count data


count_data.data = readcell(scr_para.count_data);

count_data.sample_names = count_data.data(1,:);
count_data.feature_name = count_data.data(2:end,1);

count_data.data = cell2mat(count_data.data(2:end,2:end));

% bring the metadata in the same order and then add it to object 
[~, idx, idy] = intersect(count_data.sample_names(arrayfun(@(x) isstr(x{:}),count_data.sample_names)),...
                          sample_metadata.SampleID);
count_data.sample_names =count_data.sample_names(idx);
count_data.data = count_data.data(:,idx);
count_data.metadata = sample_metadata(idy,:);

count_data.pseudo_count = count_data.data + 1;
count_data.pseudo_log2 = log2(count_data.data+1);


clear sample_metadata idx idy
%% PCA plot

mapcaplot(expression_data.data',expression_data.sample_names);

%% perform QC

disp("###----- perform QC -----###") 
disp("###----- PCA")

[coeff,score,latent,tsquared,explained] = pca(count_data.pseudo_count');
scatter(score(:,1),score(:,2) )

figure
hold on
for x = unique(count_data.metadata.(scr_para.columns_to_define_model_samples_on))'
    disp(x);
    idx = contains(count_data.metadata.(scr_para.columns_to_define_model_samples_on),x{:});
    scatter(score(idx,1),score(idx,2))
end   

title('PCA')
xlabel([num2str(1), ' component: ', num2str(explained(1))])
ylabel([num2str(2), ' component: ', num2str(explained(2))])
legend(unique(count_data.metadata.(scr_para.columns_to_define_model_samples_on))' ,'location','best')
 
saveas(gcf, [scr_para.QC_figures_path  date  '_pseudocount_PCA_.png']);

%%

disp("###----- barplot expression per sample")
bar(sum(count_data.data))
title('Number of reads per sample: ')
xlabel("samples")
ylabel("# of reads")
xticklabels(count_data.sample_names)
saveas(gcf,[scr_para.QC_figures_path  date  '_expression_read_counts_barplot.png']);

bar(sum(count_data.data == 0,1))
title('Number of 0 per sample')
xlabel("samples")
ylabel("# of zeros")
xticklabels(count_data.sample_names)
saveas(gcf, [scr_para.QC_figures_path  date  '_expression_zero_counts_barplot.png']);

disp("###----- boxplot expression per sample")

figure
boxplot(count_data.data)
saveas(gcf, [scr_para.QC_figures_path  date  '_expression_unnormalized_boxplots.png']);

figure
boxplot(expression_data.data)
saveas(gcf, [scr_para.QC_figures_path  date   '_expression_normalized_boxplots.png']);

disp("###----- boxplot expression per sample - zeros")

% without zeros
data2=expression_data.data;
data2(data2==0)=NaN;
figure
boxplot(data2)

saveas(gcf, [scr_para.QC_figures_path  date   '_expression_normalized_boxplots_nozeros.png']);

clear data2

% without zeros
data2=count_data.data;
data2(data2==0)=NaN;
figure
boxplot(data2)

saveas(gcf, [scr_para.QC_figures_path  date  '_expression_unnormalized_boxplots_nozeros.png']);

clear data2


disp("###----- density plots of normalized data per sample")

figure
hold on
fpkm = count_data.data;
fpkm(fpkm==0)=NaN; %remove zeros for densityplot
fpkm=log2(fpkm); %log2 scaling
for i=1:size(fpkm,2)
    [probability_estimate,xi] = ksdensity(fpkm(:,i));
    plot(xi,probability_estimate,':k','LineWidth',1);
end
title('Distribution per sample unnormalized raw counts')
hold off
saveas(gcf, [scr_para.QC_figures_path  date   '_density_dist_unnormalized_data.png']);

figure
hold on
fpkm = expression_data.data;
fpkm(fpkm==0)=NaN; %remove zeros for densityplot
fpkm=log2(fpkm); %log2 scaling
for i=1:size(fpkm,2)
    [probability_estimate,xi] = ksdensity(fpkm(:,i));
    plot(xi,probability_estimate,':k','LineWidth',1);
end
title('Distribution per sample fpkm')
hold off
saveas(gcf, [scr_para.QC_figures_path  date   '_density_dist_normalized_fpkmdata.png']);


%% Load medium

[NUM,TXT,RAW]=xlsread(scr_para.medium_used_file);
medium_RPMI=TXT(2:end,2);
medium_RPMI_EX=cellfun(@(x)['EX_' x],medium_RPMI,'uni',false);
T=table(medium_RPMI_EX,medium_RPMI,NUM(:,3));




%% LOAD MODEL: Recon and remove unwanted export and take up reactions


load(scr_para.model_used) 
% find the intersection of reactions in the flux file, model & the medium                            
[~,idx, idx_fluxes_in_model] = intersect(T.medium_RPMI_EX,model.rxns);                         
model.lb(idx_fluxes_in_model) = -T{idx,"Var3"};

% set the rxns lower and upper bound, rxns that we set that we do not want
% to have, are those also reasonable in my case, for my data ? 
model.ub(find(ismember(model.rxns,split(scr_para.unwanted_uptakes_export_ub, ";"))))=0; 
model.lb(find(ismember(model.rxns,split(scr_para.unwanted_uptakes_export_lb, ";"))))=0; 



%% BUILD generic CONSISTENT model - fast consistency check (fastcc)

A = fastcc_4_rfastcormics(model, 1e-4, 1);

% remove non consistent reactions from model
model=removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));
% check if the biomass reactions are still there
model.rxns(find(contains(model.rxns,'biomass')))

model_orig=model;
% check if the created model is now really consistent
A = fastcc_4_rfastcormics(model, 1e-4, 1);

clear model
%% DISCRETIZE expresssion data 

expression_data.discretized = discretize_FPKM(expression_data.data, expression_data.sample_names,1); % in cases where the left distribution is higher than the right one
%discretized = discretize_FPKM(expression_data.data, expression_data.colnames,1); %with figures, will save figures in Figures folder

figure
num_disc = hist(expression_data.discretized,3);
perc_disc = (num_disc./sum(num_disc,1)) *100;
%bar(perc_disc','stacked')
bar(num_disc','stacked')

%% BUILT CONTEXT SPECIFIC MODELS -> reconstruction using rFASTCORMICS

load(scr_para.gene_dic_file)

subSys=vertcat(model_orig.subSystems{:});
 
optional_settings.unpenalized = model_orig.rxns(ismember(subSys,scr_para.unpenalizedSystems));
optional_settings.func = {'DM_atp_c_', 'biomass_reaction'}; %, T.medium_RPMI_EX{:}}; %biomass_maintenance %-> c
optional_settings.not_medium_constrained = scr_para.not_medium_constrained;
%optional_settings.medium = T.medium_RPMI; %(add media instead)
biomass_rxn = {'biomass_reaction'} 

                        
condition_models = struct();

condition_column = scr_para.columns_to_define_model_samples_on;
% get the index of the samples in every defined group
for cond = unique(expression_data.metadata.(condition_column))'
         % transform the array, the for loop loops over the rows, so if the elements over which you want to loop over are defined in cells in one row /not column then the for loop will concat all elements instead of looping over them
        idx = contains(expression_data.metadata.(condition_column),cond)

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


md_file_name = [scr_para.save_models_to  date  "_cond_models.mat"];  % Convert datetime object to string
disp(md_file_name);
dat_file_name = [scr_para.save_models_to  date  "_input_data_obj_models.mat"];  % Convert datetime object to string
disp(dat_file_name);
save(md_file_name, 'condition_models')
save(dat_file_name, 'expression_data')
orig_file_name = [scr_para.save_models_to  date  "_orig_model.mat"];  % Convert datetime object to string
disp(orig_file_name);
save(orig_file_name, 'model_orig')

