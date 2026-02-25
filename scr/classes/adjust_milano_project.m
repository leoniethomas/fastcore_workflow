% adjust milano project file 
load("milano_project_20260224_1717.mat")
% make the discretization visualization work for the structural comparison
% script - adding metadata to the model, since for my models I have
% multiple samples going into one model and therefore multiple columns
% stored in the disc and mappedrxn, I need to know the names & identifier
% which model they belong to in order to put that on the x axis in the
% barplot 
milano_project.models.SDHB_0912.sample_metadata = table({'SDHB_0912'},{'SDHB_0912'},VariableNames=["sample", "condition"]);
milano_project.models.EV_1201.sample_metadata = table({'EV_1201'},{'EV_1201'},VariableNames=["sample", "condition"]);
milano_project.models.SDHB_1201.sample_metadata = table({'SDHB_1201'},{'SDHB_1201'},VariableNames=["sample", "condition"]);
milano_project.models.EV_0912.sample_metadata = table({'EV_0912'},{'EV_0912'},VariableNames=["sample", "condition"]);

% by default the first column in the metadata table is considered to be the
% name of the sample, then you have to give the name of the column that is
% defining the model this sample belongs to - in your case that is the same
% thing since you only have one discretization column your sample name is
% your condition - therefore we only need one row 
milano_project.models.SDHB_0912.settings.script_parameters.columns_to_define_model_samples_on = "condition";
milano_project.models.SDHB_1201.settings.script_parameters.columns_to_define_model_samples_on = "condition";
milano_project.models.EV_0912.settings.script_parameters.columns_to_define_model_samples_on = "condition";
milano_project.models.EV_1201.settings.script_parameters.columns_to_define_model_samples_on = "condition";



%% deal witht the medium 

modelNames = ["SDHB_0912", "EV_0912"];

for i = 1:numel(modelNames)
    name = modelNames{i};
    
    tmp = milano_project.models.(name).settings.medium;
    
    milano_project.models.(name).settings.medium = struct();
    milano_project.models.(name).settings.medium.medium_composition = tmp;

    milano_project.models.(name).settings.medium.manual_set_boundaries.unwanted_import = ["EX_o2s[e]", "EX_h2o2[e]", "EX_ppi[e]"];
    milano_project.models.(name).settings.medium.manual_set_boundaries.wanted_import = ["EX_o2[e]", "EX_h[e]"];
    milano_project.models.(name).settings.medium.manual_set_boundaries.wanted_export = strings(0);
    milano_project.models.(name).settings.medium.manual_set_boundaries.unwanted_export = strings(0);
end

modelNames = ["SDHB_1201", "EV_1201"];

for i = 1:numel(modelNames)
    name = modelNames{i};
    
    tmp = milano_project.models.(name).settings.medium;
    
    milano_project.models.(name).settings.medium = struct();
    milano_project.models.(name).settings.medium.medium_composition = tmp;

    milano_project.models.(name).settings.medium.manual_set_boundaries.unwanted_import = ["EX_o2s[e]", "EX_h2o2[e]", "EX_ppi[e]", "EX_cit[e]", "EX_succ[e]", "EX_glu_L[e]"];
    milano_project.models.(name).settings.medium.manual_set_boundaries.wanted_import = ["EX_o2[e]", "EX_h[e]"];
    milano_project.models.(name).settings.medium.manual_set_boundaries.wanted_export = strings(0);
    milano_project.models.(name).settings.medium.manual_set_boundaries.unwanted_export = strings(0);
end


%%

%parametersAnalysis = readtable('./scr/defaultParametersAnalysis.csv');
modelList = ["SDHB_0912","EV_0912"];
modelList = ["SDHB_1201","EV_1201"];
%analysisList = ["FVA", "FBA"];%, "sampling"];

%milano_project = singleModelAnalysis(milano_project,modelList,analysisList,parametersAnalysis);

[milano_project, analysisID] = chooseActiveAnalysisForComparison(milano_project,modelList);

referenceModel = "consistent_medium_constrained_model";
comparisonAnalysisList = ["modelStructureComparison",...
                          "modelFunctionalComparison", ...
                          "modelsComparisonSampling"];%,...
                          %"IDAREoutput"];
identifier = "";


[milano_project, comparison_name] =  modelsComparison(milano_project,modelList, analysisID,...
                                             referenceModel,comparisonAnalysisList,identifier);


%% get figures for specific subsystem/set of rxn idxs 



% get ids for glycolysis
choosen_subsystem = "Glycolysis/gluconeogenesis";
idx_subsystem_reference_model = find(choosen_subsystem == string(milano_project.models.(referenceModel).model.subSystems));



% visualize intersection on basis of the reaction presence in this
% subsystem
subsystem_feature_presence = milano_project.comparisons.(comparison_name).rxn_mapping_table{idx_subsystem_reference_model,:} ~= 0;
fig = plotFlexibleVenn(subsystem_feature_presence,...
                 milano_project.comparisons.(comparison_name).modelNames, ... 
                 "Structural model comparison: rxns presence in the " + choosen_subsystem);



% get the FBA/FVA plot for the subsystem 

get_flux_plot(milano_project,comparison_name,idx_subsystem_reference_model, "FVA",true, "threshold_flux","upper")

% visualize_flux(milano_project,comparison_name,[],{idx_subsystem_reference_model},...
%                                  ["Glycolysis"]);
% 
% visualize_fluxsum(milano_project,comparison_name,[],{idx_subsystem_reference_model},...
%                                  ["Glycolysis"]);


