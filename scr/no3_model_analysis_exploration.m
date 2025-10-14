%% analysis of the context specific models 
% what are the main questions for now ? 

% growth rate ? is it reasonable ? 


%% define VARIABLES 

clearvars -except solverOK, close all, clc % clean environment
delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names


addpath(genpath("C:\Users\leonie.thomas\plot2svg"))
addpath(genpath("C:\Users\leonie.thomas\rFASTCORMICS"))
addpath(genpath("C:\Users\leonie.thomas\scFASTCORMICS"))
changeCobraSolver("ibm_cplex");



%% define script parameters
tic 
model_id = "20250716_0612";
work_on_fastcore_exp = 0;
%project_path = "\\atlas.uni.lux\fstc_sysbio\0- UserFolders\Leonie.THOMAS\projects\20250225_glynn_bulk_metabolic_model";
project_path = "/Volumes/FSTC_SYSBIO/0- UserFolders/Leonie.THOMAS/projects/20250225_glynn_bulk_metabolic_model";
path_to_model_to_analyse = project_path + filesep + "context_specific_models" + filesep + model_id;
cd (project_path)
addpath(genpath(project_path))

input_paramters = dir(path_to_model_to_analyse + filesep + "*def_run_paramters.txt");
input_paramters = [input_paramters.folder filesep input_paramters.name];
input_paramters = readtable(input_paramters);
scr_para = cell2struct(input_paramters{:,"value"}, input_paramters{:,"slot_name"});
scr_para.model_to_load = model_id + "_cond_models.mat";
scr_para.model_workspace_to_load = model_id + "_workspace_cond_models.mat";


toc

%% load the created models with their whole workspace
% - for older versions there is not fastcore_exp file - therefore the
% conidtion and workspace file are read in 

if work_on_fastcore_exp
    load(path_to_model_to_analyse + filesep +   model_id + "_fastcore_exp.mat") % load the condition specific models created with rFASTCORMICS
else
    load(path_to_model_to_analyse + filesep +   model_id + "_cond_models.mat")
    load(path_to_model_to_analyse + filesep +   model_id + "_workspace_cond_models.mat")
    
    exp = fastcore_experiment(model_orig,dico);
    exp.condition_models = condition_models;
end

scr_para.results_path = project_path + filesep + "analysis" + filesep + model_id;
scr_para.objective = 'biomass';
scr_para.remove_unused_genes = 1;

mkdir(scr_para.results_path);
results = struct();


%writeCbModel(condition_models.MDA_MB231_Cont_NO,'format', 'json','fileName','model_Cont_NO.json')
%writeCbModel(condition_models.MDA_MB231_Cont_VC,'format', 'json','fileName','model_Cont_VC.json')

%% 
analysis_results = model_analysis(exp); 

%% plot jaccard similarity score for rxn presence 

[fig,J] = get_jaccard_similarity(analysis_results)

saveas(fig,scr_para.results_path + "\rxn_occurence_jaccard_distance.png");
results.jaccard = J;
clear J


%% visualize intersections rxn presence 

% get indices of outer and intersections
idx = analysis_results.get_intersection_plot("reaction_presence")

%% Pathway analysis

% get count of rxn per pathway per model -> pathway presence
analysis_results = analysis_results.get_pathway_counts(exp);
% visualize pathway presence
analysis_results.get_pathway_plot(1:15,"")

% get the subsystems for the outersection rxn

%% scatter plot based on pathway activity 

models_to_compare = string(analysis_results.model_names(1:2));
plotting_data = analysis_results.get_pathway_plot(1:200,"relative");
labeling_data = analysis_results.get_pathway_plot(1:20,"relative");

figure
scatter(plotting_data.(models_to_compare(1)), plotting_data.(models_to_compare(2)))
hold on 
text(labeling_data.(models_to_compare(1)), labeling_data.(models_to_compare(2)),...
                         labeling_data.Properties.RowNames, ...
                         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
                     
%% FBA 

exp.condition_models = structfun(@(x) changeObjective(x,'biomass_reaction'),...
                             condition_models,'UniformOutput',false);
exp.fba        = structfun(@(x) optimizeCbModel(x,'max','zero'),...
                                 exp.condition_models,'UniformOutput',false);
                         
fba_flux_matrix = exp.join_fba_output();

ex_rxns = exp.original_model.rxns(find(findExcRxns(exp.original_model)));
fba_exchange = fba_flux_matrix(ex_rxns(find(ismember(ex_rxns, fba_flux_matrix.Properties.RowNames))),:);
fba_exchange = fba_exchange(find(sum(abs(fba_exchange{:,:}),2) > 0),:)
fba_exchange.rxns = fba_exchange.Properties.RowNames;

writetable(fba_exchange,scr_para.results_path + filesep +  "consumption_export.xlsx") 

                         
disp("display jaccard distance based on which metabolites are ...")
disp("consumed:")
J = 1-squareform(pdist((fba_exchange{:,1:end-1} > 0)','jaccard'))
disp("and exported:")
J = 1-squareform(pdist((fba_exchange{:,1:end-1} < 0)','jaccard'))    


analysis_results.get_fba_plot(fba_flux_matrix)

%% 
[max,min] = exp.join_fva_output()

%% Similarity Based on Flux Variability Analysis

[exp.minFlux, exp.maxFlux] = structfun(@(x) fluxVariability(x),...
                               exp.condition_models,'UniformOutput',false);         
[exp.fva.minFlux,exp.fva.maxFlux] = join_fva_output(exp);


for y=1:length(fieldnames(exp.condition_models))
    for x=1:length(fieldnames(exp.condition_models)) 
        if x ~= y
            [overallSim, rxnSim] = FVAsimilarity([exp.fva.minFlux{:,y}, exp.fva.maxFlux{:,y}],...
                                                 [exp.fva.minFlux{:,x}, exp.fva.maxFlux{:,x}]);

            exp.fva_similarity_rxns{y,x} = rxnSim;
            exp.fva_similarity{y,x} = overallSim; 
        else
            exp.fva_similarity_rxns{y,x} = 1;
            exp.fva_similarity{y,x} = 1;
        end
    end
end

analysis_results.get_fva_sim_plot(exp.fva_similarity)


%% single gene deletion - essential genes - enrichment of essential genes

threshold = 0.5;

condition_models = structfun(@(x) changeObjective(x,'biomass_reaction'),...
                             condition_models,'UniformOutput',false);
                  
for x = fieldnames(condition_models)'
    
    modell = condition_models.(string(x));
    
    % perform the gene deletion to see which genes are essential 
    [modell.grRatio, modell.grRateKO, ...
     modell.grRateWT, ~,modell.essential_genes_del_Rxns, ~] = singleGeneDeletion(modell, 'FBA', string(modell.genes), 0, 0);
    modell.geneList = modell.genes;
    modell.essential_genes = modell.grRatio <= threshold;

    modell.essential_genes_Symbols = modell.geneList(modell.essential_genes); %get the identifiers for the essential genes
    [~,ia,ib] = intersect(regexprep(modell.essential_genes_Symbols,".1$",""), dico.ENTREZ); 
    modell.essential_genes_Symbols(ia) = dico.SYMBOL(ib); %extract the symbols
    
    condition_models.(x{:}) = modell;
end


% 
% % visualize enrichment ? what are those gene sets ? 
% 
figure
hold on
p = structfun(@(x) plot(sort(x.grRatio,'ascend')),condition_models)
xlabel('genes sorted ascending')
ylabel('growth Rate Ratio KO/WT')
legend(model_names)
saveas(gcf,scr_para.results_path + "\ess_genes_grRateKO_WT.png");
hold off;
% 
figure
hold on
p = structfun(@(x) plot(sort(x.grRateKO,'ascend')),condition_models)
xlabel('genes sorted ascending')
ylabel('growth Rate Ratio KO/WT')
legend(model_names)
saveas(gcf,scr_para.results_path + "\ess_genes_grRate.png");
hold off;

clear ia ib I p threshold modell
%% save gene essentiallity 

ess_genes = table(condition_models.MDA_MB231_Cont_NO.essential_genes,condition_models.MDA_MB231_Cont_VC.essential_genes,...
                  'RowNames', regexprep(string(condition_models.MDA_MB231_Cont_NO.genes),".1$",""),...
                  'VariableNames',string(fieldnames(condition_models))');

writetable(ess_genes,scr_para.results_path + "\ess_genes.csv",'WriteRowNames',true)

%% get all the essential genes form all the models

essential_genes = struct2array(structfun(@(x) x.essential_genes,condition_models,'UniformOutput',false));
essential_genes_del_Rxns = struct2array(structfun(@(x) x.essential_genes_del_Rxns,condition_models,'UniformOutput',false));

essential_genes_unique = essential_genes;
essential_genes_unique(find(sum(essential_genes,2)==2),:) = 0;
unique_essential_rxns_NO = essential_genes_del_Rxns(find(essential_genes_unique(:,1)),1);
unique_essential_rxns_NO_idx = find(ismember(model_orig.rxns,vertcat(unique_essential_rxns_NO{:})));
unique_essential_rxns_VC = essential_genes_del_Rxns(find(essential_genes_unique(:,2)),2);
unique_essential_rxns_VC_idx = find(ismember(model_orig.rxns,vertcat(unique_essential_rxns_VC{:})));


J = squareform(pdist(essential_genes','jaccard'));
%Jaccard similarity plots for sample models 7
fig = fun.plot_clustergram(1-J,...
                     model_names,...
                     model_names,...
                     {'Essential gene similarity based on Jaccard distance'},...
                     [100 100 800 600],...
                     altcolor);
saveas(fig,scr_para.results_path + "\essential_gene_similarity_jaccard_score.png");




%% visualize essential genes

id_genes_ess = find(sum(essential_genes,2) ~= 0);
essential_genes_non_zero = double(essential_genes(id_genes_ess,:));

[~,ib] = ismember(condition_models.MDA_MB231_Cont_NO.geneList(id_genes_ess)', dico.ENTREZ);

fig = fun.plot_clustergram(essential_genes_non_zero,...
                     dico.SYMBOL(ib)',...
                     model_names,...
                     {'essential genes per model'},...
                     [100 100 800 600],...
                     altcolor);
saveas(fig,scr_para.results_path + "\essential_genes.png");

clear J ib id_genes_ess essential_genes_non_zero essential_genes

%% gene set enrichment of essential genes
% get the top 15 terms and visulize them for all the models

enrichment = cell2mat(struct2array(structfun(@(x) x.enrichment.enrichment,condition_models,'UniformOutput',false)));

choose_gene_sets = find((abs(min(enrichment') - max(enrichment')) > 0.4)' | (sum(enrichment,2)>0.6));
gene_set_names = arrayfun(@(x) regexprep(x,"_","\_"),condition_models.MDA_MB231_Cont_NO.enrichment.("Database/website"));

fig = fun.plot_clustergram(enrichment(choose_gene_sets,:),...
                     gene_set_names(choose_gene_sets)',model_names,...
                     {'enrichment of essential genes per model in gene sets'},...
                     [100 100 800 600],...
                     altcolor);
saveas(fig,scr_para.results_path + "\enrichment_genesets_essential_genes.png");

clear choose_gene_sets enrichment gene_set_names
% 
% %% essential genes detected in which pathways
% 
% essential_genes_pathways = structfun(@(x) findRxnsFromGenes(x,x.genes(x.essential_genes))',condition_models,'UniformOutput',false);
% %but the pathways into one matrix for every model
% t = structfun(@(x) unique(struct2array(structfun(@(y) y(:,3)',x,'UniformOutput',false)')),essential_genes_pathways,'UniformOutput',false);
% uniq_pathways = unique(struct2array(t));
% affected_pathways_ess_genes = double(struct2array(structfun(@(x)ismember(uniq_pathways,x)',t,'UniformOutput',false)));
% 
% fig = fun.plot_clustergram(affected_pathways_ess_genes,uniq_pathways',model_names,...
%                      {'pathways for which essential genes were detected per model'},...
%                      [100 100 800 600],...
%                      altcolor);
% saveas(fig,scr_para.results_path + "\occurence_pathways_essential_genes.png");
% 
% clear uniq_pathways affected_pathways_ess_genes essential_genes_pathways t 
% 
% %% Drug deletion
% % define a list of drugs 
% load(scr_para.gene_drug_relation_file);
% DrugList = unique(GeneDrugRelations.DrugName);
% condition_models = structfun(@(x) changeObjective(x,'biomass_reaction'),...
%                              condition_models,'UniformOutput',false);
%                          
%                          
% [grRatio, grRateKO, grRateWT] = structfun(@(x) DrugDeletion(x,'FBA',DrugList),...
%                                           condition_models,'UniformOutput',false);
%                                       
% drug_deletion_res = struct("grRatio",struct2array(grRatio),...
%                            "grRateKO",struct2array(grRateKO),...
%                            "grRateWT",struct2array(grRateWT));
% 
% drugidxs_with_an_effect = find(sum(drug_deletion_res.grRatio,2)<(size(drug_deletion_res.grRatio,2) - 0.0001));
% fig = fun.plot_clustergram(double(drug_deletion_res.grRatio(drugidxs_with_an_effect,:)),...
%                      DrugList(drugidxs_with_an_effect)',model_names,...
%                      {'Drugdeletion - grRatio KO/WT'},...
%                      [100 100 800 600],...
%                      altcolor);
% saveas(fig,scr_para.results_path + "\drug_deletion_growthRateKO_WT.png");
% 
% results.drug_deletion = drug_deletion_res;
% 
% clear DrugList grRatio grRateKO grRateWT drug_deletion_res drugidxs_with_an_effect

%% how many no rxns ? in the different models 
model = condition_models.MDA_MB231_Cont_NO;
no_rxns =find(model.S(find(matches(model.mets,"no[c]")),:));
model.rxns(no_rxns)
formulas = printRxnFormula(model);
formulas(no_rxns)

phe_rxns = string(model.rxns(find(model.S(find(matches(model.mets,"pi[c]")),:))))
arg_rxns = string(model.rxns(find(model.S(find(matches(model.mets,"phe_L[c]")),:))))
cit_rxns = string(model.rxns(find(model.S(find(matches(model.mets,"citr_L[c]")),:))))

