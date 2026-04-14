%% Model comparison script 
% This script runs a set of default analysis for a choosen set of models.
% For these models a SingleModelAnalysis hast to be run beforehand! 

%% Script Setup


delete clone*.log % delet old log file 
feature astheightlimit 2000 % enable long file names
addpath(genpath("C:\Users\leonie.thomas\rFASTCORMICS")) % add rFASTCORMICS function to your path
changeCobraSolver("gurobi");

% download the git & set the folder which is gonna be your working directory
working_path = "/Users/leonie.thomas/Documents/fastcore_workflow_with_vanille";
cd (working_path)
addpath(genpath(working_path))
addpath(genpath("C:\Users\leonie.thomas\looplessFluxSampler"))
rmpath('/Users/leonie.thomas/cobratoolbox/src/analysis/thermo/thermoFBA')

% load your singleModel project object
%load(working_path + filesep + "context_specific_models" + filesep + "20260119_1042" + filesep + "project_23012026_1453_28012026_1508_obj_vanille_sampling_20260306_sampling.mat")
%load(working_path + filesep + "context_specific_models" + filesep + "20260326_0311" + filesep +  "20260326_0311_project.mat")
load(working_path + filesep + "context_specific_models" + filesep + "20260326_0311" + filesep + "20260326_0311_project_sampling_0704.mat")



%% perform new single cell analysis - if not already done for the given object
opts = detectImportOptions('./scr/defaultParametersAnalysis.csv');
opts.VariableTypes{3} = 'char'; % making sure that the last column with the values is read in as a character 
% Ensure all rows are read (no early stopping)
opts.DataLines = [1 Inf];
parametersAnalysis = readtable('./scr/defaultParametersAnalysis.csv', opts);

modelList = ["WT","KO","PLV"];
analysisList = ["FBA", "FVA", "sampling"];

project = singleModelAnalysis(project,modelList,analysisList,parametersAnalysis);

[project, analysisID] = chooseActiveAnalysisForComparison(project,modelList);

save(working_path + filesep + "context_specific_models" + filesep + "20260326_0311" + filesep + "20260326_0311_project_sampling_1404_loopless.mat",'project','-v7.3')

%% Main analysis 

referenceModel = "consistent_medium_constrained_model";
comparisonAnalysisList = ["modelStructureComparison",...
                          "modelFunctionalComparison"]%,...
                          %"modelsComparisonSampling"];%,...
                          %"IDAREoutput"];
identifier = "";

[project,comparison_name] = modelsComparison(project,modelList, analysisID,...
                                            referenceModel,comparisonAnalysisList,identifier);

%%

visualize_sampling_landscape(project,comparison_name, rxn_to_visualize,options)

%% venn for all reaactions ? 
subsystem_feature_presence = project.comparisons.(comparison_name).rxn_mapping_table ~= 0;

plotFlexibleVenn(subsystem_feature_presence{:,:},...
                 project.comparisons.(comparison_name).modelNames, ... 
                 "Structural model comparison: rxns presence");

%% 
coa = find(matches(string(project.models.(referenceModel).model.subSystems),"CoA synthesis"));
gly = find(matches(string(project.models.(referenceModel).model.subSystems),"Glycolysis/gluconeogenesis"));
oxpho = find(matches(string(project.models.(referenceModel).model.subSystems),"Oxidative phosphorylation"));
tca = find(matches(string(project.models.(referenceModel).model.subSystems),"Citric acid cycle"));
arg_pro = find(matches(string(project.models.(referenceModel).model.subSystems),"Arginine and proline metabolism"));

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{coa,gly,oxpho},...
                                 ["CoA synthesis", "Glycolysis/gluconeogenesis", "Oxidative phosphorylation"],...
                                 "heatmap",true);

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{coa,gly},...
                                 ["CoA synthesis", "Glycolysis/gluconeogenesis", "Oxidative phosphorylation"],...
                                 "heatmap_sample_all_features",true);



fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{coa,gly,oxpho},...
                                 ["CoA synthesis", "Glycolysis/gluconeogenesis", "Oxidative phosphorylation"],...
                                 "heatmap_sample_all_features",true,true,"ordered_fba", "reactions");


fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{coa,gly,oxpho},...
                                 ["CoA synthesis", "Glycolysis/gluconeogenesis", "Oxidative phosphorylation"],...
                                 "heatmap",true,true,"ordered_fba", "reactions");


%% subsystem of interest for Letizia 

gly = find(matches(string(project.models.(referenceModel).model.subSystems),"Glycolysis/gluconeogenesis"));
tca = find(matches(string(project.models.(referenceModel).model.subSystems),"Citric acid cycle"));
PPP = find(matches(string(project.models.(referenceModel).model.subSystems),"Pentose phosphate pathway"));
pyr = find(matches(string(project.models.(referenceModel).model.subSystems),"Pyruvate metabolism"));
purine = find(contains(string(project.models.(referenceModel).model.subSystems),"Purine "));

pyrimidine =find(contains(string(project.models.(referenceModel).model.subSystems),"Pyrimidine "));

nuc = find(matches(string(project.models.(referenceModel).model.subSystems),"Nucleotide interconversion"));
glut = find(matches(string(project.models.(referenceModel).model.subSystems),"Glutamate metabolism"));
Urea_cycle = find(matches(string(project.models.(referenceModel).model.subSystems),"Urea cycle"));

proline = find(matches(string(project.models.(referenceModel).model.subSystems),"Arginine and proline metabolism"));

% pick amino acids and lipids as one system
subs = string(project.models.(referenceModel).model.subSystems);
mask = contains(lower(subs), ...
    ["alanine","glycine","valine","leucine","isoleucine","serine","threonine","cysteine","methionine","aspartate","asparagine","glutamate","glutamine","arginine","proline","histidine","phenylalanine","tyrosine","tryptophan"]);
AA = find(mask);

lipid_subsystems = [
    "Fatty acid oxidation"
    "Fatty acid synthesis"
    "Glycerophospholipid metabolism"
    "Sphingolipid metabolism"
    "Cholesterol metabolism"
];

mask = ismember(subs, lipid_subsystems);
Lipids = find(mask);

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{gly,tca,PPP, pyr,purine,pyrimidine,nuc,glut,Urea_cycle,proline, AA},...
                                 ["gly","tca","PPP", "pyruvate","purine metabolism","pyrimidine metabolism", "nucleotide metabolism","glutamate metabolism","urea cycle" ,"proline and arginine metabolism", "Amino acid metabolism"],...
                                 "heatmap",true,true,"ordered_fba", "reactions");
fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{gly,tca,PPP, pyr,purine,pyrimidine,nuc,glut,Urea_cycle,proline, AA},...
                                 ["gly","tca","PPP", "pyruvate","purine metabolism","pyrimidine metabolism", "nucleotide metabolism","glutamate metabolism","urea cycle" ,"proline and arginine metabolism", "Amino acid metabolism"],...
                                 "heatmap",true,true,"ordered_fba", "incoming");
fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{gly,tca,PPP, pyr,purine,pyrimidine,nuc,glut,Urea_cycle,proline, AA},...
                                 ["gly","tca","PPP", "pyruvate","purine metabolism","pyrimidine metabolism", "nucleotide metabolism","glutamate metabolism","urea cycle" ,"proline and arginine metabolism", "Amino acid metabolism"],...
                                 "heatmap_sample_all_features",true,true,"ordered_fba", "reactions");

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{gly,tca,PPP, pyr,purine,pyrimidine,nuc,glut,Urea_cycle,proline, AA},...
                                 ["gly","tca","PPP", "pyruvate","purine metabolism","pyrimidine metabolism", "nucleotide metabolism","glutamate metabolism","urea cycle" ,"proline and arginine metabolism", "Amino acid metabolism"],...
                                 "heatmap",true,true,"ordered_samples", "reactions");
fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{gly,tca,PPP, pyr,purine,pyrimidine,nuc,glut,Urea_cycle,proline, AA},...
                                 ["gly","tca","PPP", "pyruvate","purine metabolism","pyrimidine metabolism", "nucleotide metabolism","glutamate metabolism","urea cycle" ,"proline and arginine metabolism", "Amino acid metabolism"],...
                                 "heatmap_sample",true,true,"ordered_samples", "incoming");

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{gly,tca,PPP, pyr,purine,pyrimidine,nuc,glut,Urea_cycle,proline, AA},...
                                 ["gly","tca","PPP", "pyruvate","purine metabolism","pyrimidine metabolism", "nucleotide metabolism","glutamate metabolism","urea cycle" ,"proline and arginine metabolism", "Amino acid metabolism"],...
                                 "heatmap",true);

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{gly,tca,PPP, pyr,purine,pyrimidine,nuc,glut,Urea_cycle,proline, AA},...
                                 ["gly","tca","PPP", "pyruvate","purine metabolism","pyrimidine metabolism", "nucleotide metabolism","glutamate metabolism","urea cycle" ,"proline and arginine metabolism", "Amino acid metabolism"],...
                                 "heatmap_sample",true);

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{gly,tca,PPP, pyr,purine,pyrimidine,nuc,glut,Urea_cycle,proline, AA},...
                                 ["gly","tca","PPP", "pyruvate","purine metabolism","pyrimidine metabolism", "nucleotide metabolism","glutamate metabolism","urea cycle" ,"proline and arginine metabolism", "Amino acid metabolism"],...
                                 "heatmap",true,true,"ordered_samples", "incoming");
%% Investigations up to the user (for example for specific subsystems)


% run a venn for a specific subsystem
choosen_subsystem = "Glycolysis/gluconeogenesis";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(referenceModel).model.subSystems));
% create the rxns presence table from it
subsystem_feature_presence = project.comparisons.(comparison_name).rxn_mapping_table{idx_subsystem_reference_model,:} ~= 0;
fig = plotFlexibleVenn(subsystem_feature_presence,...
                 project.comparisons.(comparison_name).modelNames, ... 
                 "Structural model comparison: rxns presence in the " + choosen_subsystem);

choosen_subsystem = "Exchange rxns";
idx_subsystem_reference_model = find(findExcRxns(project.models.(referenceModel).model));

fluxsum_sets = visualize_flux(project,comparison_name,[],{idx_subsystem_reference_model},...
                                 ["exchangers"]);

%%



%% genes of interest 
gene = "5831"; 
genes_in_model = string(project.models.(referenceModel).model.genes);
genes_of_interest =  genes_in_model(find(contains(genes_in_model , gene)));

findRxnsFromGenes(project.models.(referenceModel).model,char(genes_of_interest(1)))

%% get the rxns overview with FVA,FBA and reduced cost for a specified subsystem


idx_subsystem_reference_model = find(findExcRxns(project.models.(referenceModel).model));
get_flux_plot(project,"KO_vs_PLV_vs_WT__",idx_subsystem_reference_model, "FVA",true, "threshold_flux","upper")

visualize_flux(project,comparison_name,[],{idx_subsystem_reference_model},...
                                 ["import"]);

%%

choosen_subsystem = "Pentose Phosphate Pathway";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(referenceModel).model.subSystems));
get_flux_plot(project,comparison_name,idx_subsystem_reference_model, "FVA",false)


choosen_subsystem = "Glycolysis/gluconeogenesis";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(referenceModel).model.subSystems));
get_flux_plot(project,comparison_name,idx_subsystem_reference_model, "FVA",false)


choosen_subsystem = "Arginine and proline metabolism";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(referenceModel).model.subSystems));
get_flux_plot(project,comparson_name,idx_subsystem_reference_model, "FVA",true)


%% get rxns which show a reduced cost ~= 0 

replacement_value = "analysis.FBA.basis.reducedcost"; % get the fba solution values
ordered_reducedCost_matrix = getOrderedFeatureMatrix(project,project.comparisons.KO_vs_PLV_vs_WT__.modelNames,"rxns","orig_model",replacement_value);
reduced_cost_idx = find(sum(ordered_reducedCost_matrix,2) ~= 0);
get_flux_plot(project,"KO_vs_PLV_vs_WT__",reduced_cost_idx, ...
              "threshold_flux","all","FVA", true,'reducedCost',true ,...
              'title_plots',"Functional model comparison: all reactions with a ~= 0 reduced cost");
 
%% get rxns associated with a specific metabolite 

met_name_pattern = "^o2s[";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.mets, met_name_pattern, 'once')));
met_names = project.models.(referenceModel).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(referenceModel).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(referenceModel).model.rxns);

get_flux_plot(project,comparison_name,rxns_met_id, ...
              "threshold_flux","none","FVA", false ,'threshold_flux', "none", ...
              'title_plots',"Functional model comparison: all reactions including proline");

%% get rxns associated with a specific metabolite 

rxn_name_pattern = "^P5CR.*";
idx_rxn_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.rxns, rxn_name_pattern, 'once')));
rxn_names = project.models.(referenceModel).model.rxns(idx_rxn_matches);

[~,rxns_id] = ismember(rxn_names, project.models.(referenceModel).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_id, ...
              "threshold_flux","none","FVA", true ,'threshold_flux', "all", ...
              'title_plots',"Functional model comparison: all reactions including proline");

%% get rxns associated with a specific metabolite 

met_name_pattern = "^pro_L[.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.mets, met_name_pattern, 'once')));
met_names = project.models.(referenceModel).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(referenceModel).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(referenceModel).model.rxns);

get_flux_plot(project,comparison_name,rxns_met_id, ...
              "threshold_flux","none","FVA", false ,...
              'title_plots',"Functional model comparison: all reactions including pyruvate");

%% get rxns associated with a specific metabolite 

met_name_pattern = "^akg[.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.mets, met_name_pattern, 'once')));
met_names = project.models.(referenceModel).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(referenceModel).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(referenceModel).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","all","FVA", true ,...
              'title_plots',"Functional model comparison: all reactions including succinate");

%% get rxns associated with a specific metabolite 

met_name_pattern = "^galgluside_hs.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.mets, met_name_pattern, 'once')));
met_names = project.models.(referenceModel).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(referenceModel).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(referenceModel).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","none","FVA", true ,...
              'title_plots',"Functional model comparison: all reactions including lactate");

%% get rxns associated with a specific metabolite 

met_name_pattern = "^pro_L.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.mets, met_name_pattern, 'once')));
met_names = project.models.(referenceModel).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(referenceModel).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(referenceModel).model.rxns);

get_flux_plot(project,"WT_vs_KO_vs_PLV__",rxns_met_id, ...
              "threshold_flux","none","FVA", false ,...
              'title_plots',"Functional model comparison: all reactions including lactate");

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{rxns_met_id},...
                                 ["lactate"],"violin",true,false);
fluxsum_sets = visualize_flux(project,comparison_name,[],{rxns_met_id},...
                                 ["lactate"]);


%%

rxn_color = ["EX_glc[e]", "EX_gln_L[e]", "EX_pyr[e]","EX_so3[e]", "EX_retfa[e]"];

visualize_sampling_landscape(project,comparison_name,[1 2],rxn_color)

%% Glycolysis specific glucose and pyruvate


coa = find(matches(string(project.models.(referenceModel).model.subSystems),"CoA synthesis"));
gly = find(matches(string(project.models.(referenceModel).model.subSystems),"Glycolysis/gluconeogenesis"));
oxpho = find(matches(string(project.models.(referenceModel).model.subSystems),"Oxidative phosphorylation"));
tca = find(matches(string(project.models.(referenceModel).model.subSystems),"Citric acid cycle"));
arg_pro = find(matches(string(project.models.(referenceModel).model.subSystems),"Arginine and proline metabolism"));

fluxsum_sets = visualize_flux(project,comparison_name,[],{arg_pro},...
                                 ["Arginine and proline metabolism"]);

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{arg_pro},...
                                 ["Arginine and proline metabolism"],"violin",false,false);

%%

[exchange,uptake] = find(findExcRxns(project.models.(referenceModel).model));
fluxsum_sets = visualize_flux(project,comparison_name,[],{idx_subsystem_reference_model},...
                                 ["uptake"]);




%% 

met_name_pattern = "^akg.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.mets, met_name_pattern, 'once')));
met_names = project.models.(referenceModel).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(referenceModel).model, met_names);
[~,akg] = ismember(rxns_met, project.models.(referenceModel).model.rxns);

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{akg},...
                                 ["lactate"],"violin",true,true);

fluxsum_sets = visualize_flux(project,comparison_name,[],{akg},...
                                 ["lactate"]);


%%
met_name_pattern = "^gluside_hs.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.mets, met_name_pattern, 'once')));
met_names = project.models.(referenceModel).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(referenceModel).model, met_names);
[~,pro] = ismember(rxns_met, project.models.(referenceModel).model.rxns);

fluxsum_sets = visualize_flux(project,comparison_name,[],{pro},...
                                 ["lactate"]);
%% pathway fluxsum with 

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{coa,gly,oxpho},...
                                 ["CoA synthesis", "Glycolysis/gluconeogenesis", "Oxidative phosphorylation"],...
                                 "heatmap",true);

% or without the coenzymes

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{coa,gly,oxpho},...
                                 ["CoA synthesis", "Glycolysis/gluconeogenesis", "Oxidative phosphorylation"],...
                                 "heatmap",false);

%% pathway fluxsum with predefined metabolites!

pathways_with_metabolites = get_essential_pathway_metabolites(project,"orig_model");

pathways_idx = structfun(@(x)find(matches(project.models.(referenceModel).model.mets,x)),...
                         pathways_with_metabolites, 'UniformOutput', false);
fields = fieldnames(pathways_idx);
idx_cell = cellfun(@(f) pathways_idx.(f), fields, 'UniformOutput', false);    
fluxsum_sets = visualize_fluxsum(project,comparison_name,[],idx_cell,...
                                 fields);



%% summary figure letizia 

% increase in proline export ? 
% - how does proline get transported out of the cell -> do we have a EX_pro ? 
% -> 'EX_pro_L[e]' is thrown out of the model for some reason 
% -> instead the model uses 4hpro_LT -> pro + aKG -> succ + 4hproLT
% -> put the wrong direction 

% decrease in lactate excretion with KO  ? 

% increase of glucose consumption ? 

% decrease in NADP+ flux through the PPP

% decrease in the activitiy of arginine to citrulline 


% nucleotide biosynthesis goes down 

% AKGDm -> decrease in activity-> upon KO -> According to FVA and FBA -> 


%% arginine and proline metabolism

choosen_subsystem = "Arginine and proline metabolism";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(referenceModel).model.subSystems));
get_flux_plot(project,comparison_name,idx_subsystem_reference_model, "FVA",true)

ppp = find(matches(string(project.models.(referenceModel).model.subSystems),choosen_subsystem));

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{ppp},...
                                 [choosen_subsystem]);

fluxsum_sets = visualize_flux(project,comparison_name,[],{ppp},...
                                 [choosen_subsystem]);
%% pentose phosphate pathway 

choosen_subsystem = "Glycolysis/gluconeogenesis";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(referenceModel).model.subSystems));


ppp = find(matches(string(project.models.(referenceModel).model.subSystems),choosen_subsystem));

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{ppp},...
                                 [choosen_subsystem]);

fluxsum_sets = visualize_flux(project,comparison_name,[],{ppp},...
                                 [choosen_subsystem]);

%% pentose phosphate pathway 

choosen_subsystem = "Pentose phosphate pathway";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(referenceModel).model.subSystems));
get_flux_plot(project,"KO_vs_WT__",idx_subsystem_reference_model, "FVA",false)

ppp = find(matches(string(project.models.(referenceModel).model.subSystems),choosen_subsystem));

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{ppp},...
                                 [choosen_subsystem]);

fluxsum_sets = visualize_flux(project,comparison_name,[],{ppp},...
                                 [choosen_subsystem]);

%% TCA and alpha ketoglutarate 


choosen_subsystem = "Citric acid cycle";
% pull the subsystem presence from the stored rxns mapping table
idx_subsystem_reference_model = find(choosen_subsystem == string(project.models.(referenceModel).model.subSystems));
get_flux_plot(project,"KO_vs_PLV_vs_WT__",idx_subsystem_reference_model, "FVA",true)

tca = find(matches(string(project.models.(referenceModel).model.subSystems),choosen_subsystem));

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{tca},...
                                 [choosen_subsystem]);

fluxsum_sets = visualize_flux(project,comparison_name,[],{tca},...
                                 [choosen_subsystem]);


%% LACTATE FIGURE 


% FBA values lactate 

met_name_pattern = "^lac_L.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.mets, met_name_pattern, 'once')));
met_names = project.models.(referenceModel).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(referenceModel).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(referenceModel).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","all","FVA", false ,...
              'title_plots',"Functional model comparison: all reactions including lactate");


get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","none","FVA", true ,...
              'title_plots',"Functional model comparison: all reactions including lactate");

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{rxns_met_id},...
                                 ["lactate"],"violin",true,false);

fluxsum_sets = visualize_flux(project,comparison_name,[],{rxns_met_id},...
                                 ["lactate"]);

% FBA values pyruvate

met_name_pattern = "^pro_L[.*";
idx_met_matches = find(~cellfun(@isempty, regexp(project.models.(referenceModel).model.mets, met_name_pattern, 'once')));
met_names = project.models.(referenceModel).model.mets(idx_met_matches);

rxns_met = findRxnsFromMets(project.models.(referenceModel).model, met_names);
[~,rxns_met_id] = ismember(rxns_met, project.models.(referenceModel).model.rxns);

get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","none","FVA", true ,...
              'title_plots',"Functional model comparison: all reactions including proline");


get_flux_plot(project,"KO_vs_PLV_vs_WT__",rxns_met_id, ...
              "threshold_flux","none","FVA", true ,...
              'title_plots',"Functional model comparison: all reactions including proline");

fluxsum_sets = visualize_fluxsum(project,comparison_name,[],{rxns_met_id},...
                                 ["proline"],"violin",true,false);

fluxsum_sets = visualize_flux(project,comparison_name,[],{rxns_met_id},...
                                 ["proline"]);


%% visualize the expression values + their discretization status

exprCells = struct2cell(structfun(@(x)x.sample_metadata.condition,rmfield(project.models,setdiff(fieldnames(project.models),modelList)),'UniformOutput',false));
condition = string(cat(1, exprCells{:})); 

gene_names = string(project.models.KO.discretized_data.gene_names);
exprCells = struct2cell(structfun(@(x)x.expression_data,rmfield(project.models,setdiff(fieldnames(project.models),modelList)),'UniformOutput',false));
all_expr = cat(2, exprCells{:});  % Result: 24720 x total_samples


all_expr(find(contains( gene_names, "PDHX")),:)

exprCells = struct2cell(structfun(@(x)x.discretized_data.values,rmfield(project.models,setdiff(fieldnames(project.models),modelList)),'UniformOutput',false));
all_disc = cat(2, exprCells{:});  % Result: 24720 x total_samples


all_disc(find(contains( gene_names, "PDHX")),:)


PDHX_expr = all_expr(contains(gene_names, "PDHX"), :);

% Make boxplot grouped by condition
figure;
boxplot(PDHX_expr, condition);
ylabel('Expression of PDHX');
title('PDHX Expression Across Conditions [FPKM]');
grid on;



%% check how different fba solutions are 
modelNames = ["KO","WT", "PLV"];
reference = "orig_model";
ordered_FBA_matrix = getOrderedFeatureMatrix(project,modelNames,...
                                                             "rxns",reference,"analysis.FBA.v");

length(find(sum(ordered_FBA_matrix ~= 0,2) > 0))
length(find(sum(ordered_FBA_matrix ~= 0,2) == 1))
length(find(sum(ordered_FBA_matrix ~= 0,2) == 2))



%% ----------------------------------
%%%%%%%%%% Correlation analysis 
%%% ---------------------------------


data = project.models.WT.analysis.analysis_20260211_1507.sampling.samples;

biomass_id = find(matches(project.models.WT.model.rxns, "biomass_reaction"));
export = find(findExcRxns(project.models.WT.model) & project.models.WT.analysis.FBA.v > 0);
import = find(findExcRxns(project.models.WT.model) & project.models.WT.analysis.FBA.v < 0);
% get glucose, lactose, h2o, o2



% Inputs
biomass_flux = data(biomass_id, :);   % 1 x nSamples
rxn_fluxes   = data(import, :);       % nMet x nSamples
nrxns = size(rxn_fluxes, 1);
rxns_names = project.models.WT.model.rxns;

% Determine a roughly square layout
nCols = ceil(sqrt(nrxns));
nRows = ceil(nrxns / nCols);

figure
tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact')

for m = 1:nrxns
    nexttile
    scatter(biomass_flux, rxn_fluxes(m,:), 15, 'filled', 'MarkerFaceAlpha',0.6)
    xlabel('Biomass flux')
    ylabel('Metabolite flux')
    if exist('met_names','var')
        title(rxns_names(import(m)), 'Interpreter','none', 'FontSize',9)
    else
        title(sprintf('Metabolite %d', m), 'FontSize',9)
    end
    grid on
end

% Optional: link axes for all tiles
ax = findall(gcf,'Type','axes');
linkaxes(ax, 'x')  % link x-axis so all tiles share biomass scale

% Optional overall title
sgtitle('Metabolite fluxes vs biomass')






%% correlation analysis
% check first which rxns are structurally coupled
[reduced_net, fctable, blocked] = QFCA(project.models.WT.model,0, "gurobi");
%% just for one model at a time
% then perform correlation analysis between for the sampling
[correlation_coeff_matrix,pval] = corr(project.models.WT.analysis.sampling.samples','Type','Spearman');

% correlation high coefficient, highly significant + check that they are
% structurally coupled
correlated_coulpled_rxns = abs(correlation_coeff_matrix) > 0.9 & pval < 0.001 & fctable ~= 0;

%% which rxns are correlated + coupled to the biomass reaction 
biomass_id = find(matches(project.models.WT.model.rxns, "biomass_reaction"));

rxn_id = find(correlated_coulpled_rxns(biomass_id,:));

% Inputs
data = project.models.WT.analysis.sampling.samples;   % nRxns x nSamples
biomass_flux = data(biomass_id, :);
rxn_names = project.models.WT.model.rxns(rxn_id);

nRxn = numel(rxn_id);

% Choose a roughly square layout
nCols = ceil(sqrt(nRxn));
nRows = ceil(nRxn / nCols);

figure
tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact')

for k = 1:nRxn
    rxn_name = regexprep(rxn_names{k},"_"," ");
    nexttile
    scatter( ...
        data(rxn_id(k), :), ...
        biomass_flux, ...
        12, ...
        'filled', ...
        'MarkerFaceAlpha', 0.6 ...
    );

    ylabel('Biomass flux')
    xlabel(['Reaction flux '  rxn_name])

    title(rxn_name, 'Interpreter','none', 'FontSize', 9)
    grid on
end

% Optional overall title
sgtitle('Correlated reactions vs biomass (WT sampling)')



%% adjacency - shortest path

% INPUTS 
% S              : stoichiometric matrix (nMets x nRxns)
% corrMat        : reaction–reaction correlation matrix from sampling (nRxns x nRxns)
% rxn_idx_corr   : indices of reactions you identified as highly correlated
% rxn_count      : number of reactions each metabolite participates in (for currency detection)

% STEP 1 — Remove currency / coenzyme metabolites
% Highly connected metabolites (ATP, NAD, CoA, etc.) create artificial shortcuts
% Threshold must be tuned per model (50–100 is typical for Recon-scale models)

currency_threshold = 50;
currency_mets = rxn_count > currency_threshold;

% Build binary metabolite–reaction participation matrix
B = S ~= 0;

% Remove currency metabolites BEFORE building the graph
B(currency_mets, :) = false;

% STEP 2 — Build reaction–reaction adjacency matrix
% Two reactions are connected if they share at least one (non-currency) metabolite

A = (B' * B) > 0;     % rxns x rxns adjacency
A(eye(size(A))==1) = 0;  % remove self-loops

% Create graph object
G = graph(A);

% STEP 3 — Compute shortest paths between correlated reactions
% Shortest path = minimum number of reaction–reaction hops
% This captures topological proximity in the metabolic network

% Pairwise distances between all correlated reactions
D = distances(G, rxn_idx_corr, rxn_idx_corr);

% Example: shortest path between two specific correlated reactions
i = rxn_idx_corr(1);
j = rxn_idx_corr(2);

[path, dist] = shortestpath(G, i, j);
% path = sequence of reaction indices connecting i → j
% dist = path length (number of edges)

% STEP 4 — Interpret distances (conceptual guide)
% dist = 1–3  → local pathway coupling
% dist = 4–6  → mid-range network coordination
% dist = Inf  → correlation not explainable by network topology
%
% IMPORTANT:
% - Short paths do NOT prove mechanistic coupling
% - They provide topological context for correlations
% - Always interpret together with qFCA results

% STEP 5 — Optional visualization of a shortest path
subG = subgraph(G, path);
figure
plot(subG, 'Layout','layered')
title('Shortest metabolic path between correlated reactions')

% KEY TAKEAWAY (conceptual, not code)
% qFCA  → tells you which correlations are structurally meaningful
% Sampling + correlation → tells you which relationships are behaviorally expressed
% Shortest paths → tell you whether correlations are topologically local or global
%
% Used together, these three layers give structure + behavior + network context



%% export to IDARE and escher!! 


%% export to R
