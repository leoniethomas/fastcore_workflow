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

writeCbModel(condition_models.MDA_MB231_Cont_NO,'format', 'json','fileName','model_Cont_NO.json')
writeCbModel(condition_models.MDA_MB231_Cont_VC,'format', 'json','fileName','model_Cont_VC.json')

%% script initialization

changeCobraSolver("ibm_cplex");
mkdir(scr_para.results_path);
results = struct();

% load functions
fun = functions_no3;

%% get model size - HOW BIG ARE OUR MODELS ? 

condition_models = structfun(@(x) removeUnusedGenesFastbox(x,1), ... % remove unused genes to get model size in terms of genes active 
                             condition_models,'UniformOutput',false);
                         


% simple quantities per model - #rxns #genes #metabolites in the models
results.model_size = array2table(struct2array(structfun(@(x) {numel(x.rxns);numel(x.mets);numel(x.genes)}, ...
                                                        condition_models,'UniformOutput',false))',...
                                 'VariableNames',{'count_reactions','count_metabolites','count_genes'},...
                                 'RowNames',model_names)
results.model_size


load(path_to_model_to_analyse + "\" + scr_para.model_to_load) % load context specific models again, to have all the genes in them, needed for later comparison
condition_models = rmfield(condition_models,'MDA_MB231_HERVK_C_NO')
condition_models = rmfield(condition_models,'MDA_MB231_HERVK_C_VC')
condition_models = rmfield(condition_models,'MDA_MB231_HERVK_D_NO')
condition_models = rmfield(condition_models,'MDA_MB231_HERVK_D_VC')


%% Rxn occurence similarity - HOW SIMILAR IS THE CONTENT OF THE MODELS AT HAND ? 

% get the rxns still in the modell into a matrix where the rxns which are kept are set to 1, instead of a array which entails the indices of all the rxns kept  
AA_keep = struct2array(structfun(@(x) fun.reverse_find_from_indices(length(model_orig.rxns),x.AA), ...
                       condition_models,'UniformOutput',false));
% in that format we can compute a distance matrix
J = squareform(pdist(AA_keep','jaccard'));
fig = fun.plot_clustergram(1-J,...
                     model_names,...
                     model_names,...
                     {'Model similarity based on Jaccard distance of rxns existence in the model!'},...
                     [100 100 800 600],...
                     altcolor);
saveas(fig,scr_para.results_path + "\rxn_occurence_jaccard_distance.png");
results.jaccard = J;

clear J
%% inner and outersection 
inter = length(intersect(condition_models.MDA_MB231_Cont_VC.rxns, condition_models.MDA_MB231_Cont_NO.rxns));
out_VC = length(setdiff(condition_models.MDA_MB231_Cont_VC.rxns, condition_models.MDA_MB231_Cont_NO.rxns));
out_NO = length(setdiff(condition_models.MDA_MB231_Cont_NO.rxns, condition_models.MDA_MB231_Cont_VC.rxns))


%% Pathway analysis 

a = arrayfun(@(x) fun.get_pathway_counts(AA_keep(:,x),model_orig),...
             1:size(AA_keep,2),...
             'UniformOutput', false);
        
a = cat(2,a{:});
b = fun.get_pathway_counts(ones(1,size(model_orig.subSystems,1))',model_orig);
PathwayActivity_A_keep = a(:,:)./b;
%filter_0_idx = find(sum(PathwayActivity_A_keep,2) > size(PathwayActivity_A_keep,2)*0.7);
filter_0_idx  = find(var(PathwayActivity_A_keep')' > 0);
pathways_to_plot = unique([model_orig.subSystems{:}]);
pathways_to_plot = pathways_to_plot(filter_0_idx);
b = b(filter_0_idx);
pathways_to_plot = pathways_to_plot + "  -  " + string(b)';


fig = fun.plot_clustergram(PathwayActivity_A_keep(filter_0_idx,:),...
                     pathways_to_plot,...
                     model_names,...
                     {'Pathway activity for all models [count rxn per pathway/ " for consistent model]'},...
                     [100 100 800 600],...fig
                     altcolor)
saveas(fig,scr_para.results_path + "\pathway_activity_jaccard_score.png");
plot2svg(['./analysis/' char(model_id) '/pathway_activity_jaccard_score_svg.svg'],fig)
results.pathway_activity = PathwayActivity_A_keep;


%clear AA_keep PathwayActivity_A_keep a b filter_0_idx pathways_to_plot

% get colorbar for set size 
b = [90 178 1000 40 98 31 47 93 165 211 450 133 34 40 46 8 45 32 151 71 32 65 105 27 961 15 41 14 11 16 21 164 59 19 3 7];
heatmap_data = reshape(b, [], 1); % column vector

% Create the heatmap using imagesc
imagesc(heatmap_data);


%% Pathway analysis of outersection
AA_keep_outer = AA_keep; 
AA_keep_outer(find(sum(AA_keep_outer,2) ~= 1),:) = 0 ;

a = arrayfun(@(x) fun.get_pathway_counts(AA_keep_outer(:,x),model_orig),...
             1:size(AA_keep_outer,2),...
             'UniformOutput', false);
        
a = cat(2,a{:});
b = fun.get_pathway_counts(ones(1,size(model_orig.subSystems,1))',model_orig);
%PathwayActivity_A_keep = a(:,:)./b;
PathwayActivity_A_keep = a;
%filter_0_idx = find(sum(PathwayActivity_A_keep,2) > size(PathwayActivity_A_keep,2)*0.7);
%filter_0_idx  = find(var(PathwayActivity_A_keep')' > 0);
filter_0_idx = sum(PathwayActivity_A_keep,2)>0;
pathways_to_plot = unique([model_orig.subSystems{:}]);
pathways_to_plot = pathways_to_plot(filter_0_idx);
b = b(filter_0_idx);
pathways_to_plot = pathways_to_plot + "  -  " + string(b)';


fig = fun.plot_clustergram(PathwayActivity_A_keep(filter_0_idx,:),...
                     pathways_to_plot,...
                     model_names,...
                     {'Pathway activity for all models [count rxn per pathway/ " for consistent model]'},...
                     [100 100 800 600],...fig
                     altcolor)
saveas(fig,scr_para.results_path + "\pathway_activity_outer_jaccard_score.png");
plot2svg(['./analysis/' char(model_id) '/pathway_activity_outer_jaccard_score_svg.svg'],fig)
results.pathway_activity_outer = PathwayActivity_A_keep;

%clear AA_keep PathwayActivity_A_keep a b filter_0_idx pathways_to_plot

model_orig.rxns(AA_keep_outer(:,2))

%% FBA 

condition_models = structfun(@(x) changeObjective(x,'biomass_reaction'),...
                             condition_models,'UniformOutput',false);
AA_mets_imp = zeros(numel(model_orig.rxns),length(condition_models));
AA_mets_exp = zeros(numel(model_orig.rxns),length(condition_models));

filename = 'consumption_export.xlsx';

for x=fieldnames(condition_models)'
    
    model = condition_models.(string(x));
    model.solution = optimizeCbModel(model,'max','zero');
    ex_rxns = find(findExcRxns(model));
    exc_reactions = table(ex_rxns,printRxnFormula(model, model.rxns(ex_rxns)), ...
                          model.rxns(ex_rxns), model.lb(ex_rxns),model.ub(ex_rxns),...
                          model.solution.v(ex_rxns),...
                          'VariableNames',{'rxn_id' 'met_exchanged' 'ex_rxns' 'lb' 'ub' 'flux'});
    
    disp("----------" + string(x) + " ---- Ex reactions flux < 0 --------------------")
    exc_reactions(exc_reactions.flux < 0,:)
    AA_mets_imp(:,find(matches(fieldnames(condition_models),x))) = matches(model_orig.rxns,string(exc_reactions{exc_reactions.flux < 0,'ex_rxns'}));
    exc_reactions(exc_reactions.flux < 0,:),scr_para.results_path + filesep +  filename,'Sheet',string(x) + "_consumption"
    writetable(exc_reactions(exc_reactions.flux < 0,:),scr_para.results_path + filesep +  filename,'Sheet',regexprep(string(x), "HERVK", "HERV") + "_consumption") % writetable has a problem with the K for some reason...dont know why...
    disp("----------" + string(x) + " ---- Ex reactions flux > 0 --------------------")
    exc_reactions(exc_reactions.flux > 0,:)
    writetable(exc_reactions(exc_reactions.flux > 0,:),scr_para.results_path + filesep + filename,'Sheet',string(x) + "_export")
    condition_models.(string(x)) = model;
    AA_mets_exp(:,find(matches(fieldnames(condition_models),x))) = matches(model_orig.rxns,string(exc_reactions{exc_reactions.flux > 0,'ex_rxns'}));
    
end
clear exc_reactions x model

results.FBA = [fieldnames(condition_models),struct2cell(structfun(@(x) x.solution.f,...
                             condition_models,'UniformOutput',false))]
                         
disp("display jaccard distance based on which metabolites are ...")
disp("consumed:")
J = 1-squareform(pdist(AA_mets_imp','jaccard'))
disp("and exported:")
J = 1-squareform(pdist(AA_mets_exp','jaccard'))    
                         
%% visualize FBA results 

names = string(fieldnames(condition_models));
fluxes = full(sparse(cell2mat(arrayfun(@(x) condition_models.(names(x)).AA', 1:length(names), 'UniformOutput', false))', ...
                     cell2mat(arrayfun(@(x) ones(1, numel(condition_models.(names(x)).AA)) * x, 1:length(names), 'UniformOutput', false))',...
                     cell2mat(arrayfun(@(x) condition_models.(names(x)).solution.v', 1:length(names), 'UniformOutput', false))', ...
                     length(model_orig.rxns), length(fieldnames(condition_models))));

% jaccard similarity for FBA defined fluxes

J = squareform(pdist(fluxes','jaccard'));
fig = fun.plot_clustergram(log(1-J),...
                     model_names,...
                     model_names,...
                     {'Similarity of optimal fluxes obtained via FBA [log(jaccard  similarity score)]'},...
                     [100 100 800 600],...
                     altcolor);

% [coeff,score,latent,tsquared,explained] = pca(fluxes');
% 
% figure
% hold on
% for x = names'
%     disp(x);
%     idx = find(contains(names,x));
%     disp(idx);
%     scatter(score(idx,1),score(idx,2))
% end   
% 
% title('PCA - fluxes FBA')
% xlabel([num2str(1), ' component: ', num2str(explained(1))])
% ylabel([num2str(2), ' component: ', num2str(explained(2))])
% legend(model_names ,'location','best')
% hold off

%% Similarity Based on Flux Variability Analysis

[minFlux, maxFlux] = structfun(@(x) fluxVariability(x),...
                               condition_models,'UniformOutput',false); 
            
FVA_res.minFlux = arrayfun(@(x) fun.expand_values(condition_models.(x{:}).AA,...
                                                  minFlux.(x{:}),...
                                                  length(model_orig.rxns)),...
                           fieldnames(condition_models), 'UniformOutput', false);
                       
FVA_res.minFlux = horzcat(FVA_res.minFlux{:});   

FVA_res.maxFlux = arrayfun(@(x) fun.expand_values(condition_models.(x{:}).AA,...
                                                  maxFlux.(x{:}), ...
                                                  length(model_orig.rxns)),...
                           fieldnames(condition_models), 'UniformOutput', false);
                       
FVA_res.maxFlux = horzcat(FVA_res.maxFlux{:});  

% compute similarity on basis of the FVA
FVA_sim_reactions = cell(length(fieldnames(condition_models)),length(fieldnames(condition_models)));
FVA_sim_overall = FVA_sim_reactions;

for y=1:length(fieldnames(condition_models))
    for x=1:length(fieldnames(condition_models))                       
        [overallSim, rxnSim] = FVAsimilarity([FVA_res.minFlux(:,y), FVA_res.maxFlux(:,y)],...
                                             [FVA_res.minFlux(:,x), FVA_res.maxFlux(:,x)]);
        
        FVA_sim_reactions{y,x} = rxnSim;
        FVA_sim_overall{y,x} = overallSim;      
    end
end

fig = fun.plot_clustergram(cell2mat(FVA_sim_overall),...
                     model_names,model_names,...
                     {'FVA overall Similarity'},...
                     [100 100 800 600],...
                     altcolor);
saveas(fig,scr_para.results_path + "\FVA_similarity_score.png");

results.FVA = FVA_res;
results.FVA.sim_rxns = FVA_sim_reactions;
results.FVA.sim_overall = FVA_sim_overall;

clear FVA_res FVA_sim_reactions FVA_sim_overall minFlux maxFlux overallSim rxnSim
                 
%%
FVAsimilarity([FVA_res.minFlux(:,y), FVA_res.maxFlux(:,y)],...
                                             [FVA_res.minFlux(:,x), FVA_res.maxFlux(:,x)]);
        
 
%% save results 

save(scr_para.results_path + filesep + char(datetime('now', 'Format', 'yyyyMMdd_hhss')) + "_results_analysis.mat", 'results');

  


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
    [~,ia,ib] = intersect(modell.essential_genes_Symbols, dico.ENTREZ); 
    modell.essential_genes_Symbols(ia) = dico.SYMBOL(ib); %extract the symbols

    [modell.enrichment] = GeneEnrichments(modell.essential_genes_Symbols);
    [~,I] = sort(cell2mat(modell.enrichment.enrichment));
    modell.enrichment(I,:);
        
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



