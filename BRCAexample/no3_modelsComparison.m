%% Main script n°3: Models Comparison 
% This script is showing how to compare models stored inside a project.
% Models can be compared on structural, functional and sampling aspects.
% For both functional and sampling comparison, functional tests (as FBA, FBA or sampling)
% must have been previously performed (see tutorial script n°2 for single model analysis) 
% and are required to be store inside the project in the proper format.
% 
%% INITIALIZING THE ENVIRONNEMENT
initCobraToolbox();
changeCobraSolver('gurobi');
feature astheightlimit 2000;

%% LOADING PROJECT AND PARAMETERS FOR ANALYSIS
load('20262608_BRCAProject.mat'); % available on the zenodo link
defaultParametersAnalysis = readInParamTable('data/defaultParametersTable.csv');

%% COMPUTING THE COMPARISON
modelsToCompare = {"Control", "StageI", "StageII", "StageIII"};
[BRCAProject, analysisIDs] = chooseActiveAnalysis(BRCAProject, modelsToCompare); % by default, the most recent analyses will be chosen

% choose a specific analysis
% analysisToChoose = {'analysis_20260812_2130', 'analysis_20260812_2112', 'analysis_20260825_0901', 'analysis_20260825_1857'};
% [BRCAProject, analysisIDs] = chooseActiveAnalysis(BRCAProject, modelsToCompare, analysisToChoose);

%% VISUALIZING AUTOMATICALLY GENERATED FIGURES

referenceModel = "consistentMediumConstrainedModel"; 
comparisonList = ["structuralComparison", "functionalComparison", "samplingComparison"];
compID = "test_tuto3"; 

[BRCAProject, comparisonName] = modelsComparison(BRCAProject, modelsToCompare, referenceModel, compID, comparisonList);

%% Showing the Default figures generated during comparison 

load('20262608_bigComparisonStructComparison.mat');

compName = "Control_vs_StageI_vs_StageII_vs_StageIII__full_test_finalObject";

% --- Structural Analysis Plots

% how similar are the models to one another 
% jaccard similarity over the model structure, meaning are rxns included in
% the model
showFigure(BRCAProject.comparisons.(compName).structuralComparison.plots.jaccardDist.rxns)
% most reactions are shared between the models, which is to be expected 
% what can also be observed is that the cancer models are more similar than
% to the normal model

% lets check out the size and overlap of the models, in terms of genes,
% metabolites and reactions
showFigure(BRCAProject.comparisons.(compName).structuralComparison.plots.intersections.genes)
showFigure(BRCAProject.comparisons.(compName).structuralComparison.plots.intersections.mets)
showFigure(BRCAProject.comparisons.(compName).structuralComparison.plots.intersections.rxns)
% in these venn diagramms we see in absolute numbers what is indicated by
% the jaccard similarity heatmap

% how were the genes in the model discretized ? 
showFigure(BRCAProject.comparisons.(compName).structuralComparison.plots.dataDiscretization)
% this is more a QC measure, can be observed that roughly 1/3 of all genes
% which are in the model where discretized to be active
% this gives us an indication on how many of the genes,rxns are actually
% backed up by expression data

% Next lets answer two questions: 
% - how many of the reactions in the model are in the core (along the same
% lines as the discretization_ 
% - how many of the reactions which were in the core made it into the model
% (theoretically we would like 100 percent) but in practice we are around
% 70/80 percent, this can be adjusted by setting differen thresholds 
showFigure(BRCAProject.comparisons.(compName).structuralComparison.plots.coreReactions)
% the next question then is how many of those core reactions are specific
% or are the all shared ? 
% The models share a large amount of core reactions, which makes sense but we also see that the control samples have the highest number of specific rxns, 
% while the overlab between the two cancer models are larger than between
% the cancer models and control, stage I has the least number of specific
% reactions

% where do those core reactions come from, pathway wise 
showFigure(BRCAProject.comparisons.(compName).structuralComparison.plots.coreReactionsIntersections)
% the important observation here is that the differences in core reaction
% is not only due to Transport or exchangers, but there is also other rxns
% that make up the model specific core reactions



%the next question is where do the differences between the models come from
%in terms of rxns presecence, so structural difference
showFigure(BRCAProject.comparisons.(compName).structuralComparison.plots.reactionPathwayPresence)


%%

% Is there a difference in growth between the models ? 

showFigure(BRCAProject.comparisons.(compName).functionalComparison.plots.objValue)

% Is there a difference in what the models consume and how much ? 
showFigure(BRCAProject.comparisons.(compName).functionalComparison.plots.import)

% Is there a difference in what the models export and how much ? 
showFigure(BRCAProject.comparisons.(compName).functionalComparison.plots.export)

% How similar are the models in terms of FVA boundaries ? 
% TODO: go through the functions, check what I did, are those plots usefull ? 
showFigure(BRCAProject.comparisons.(compName).functionalComparison.plots.fvaSim.overall)

showFigure(BRCAProject.comparisons.(compName).functionalComparison.plots.fvaSim.hist)



% How much is a given pathway used in our FBA 

% used in terms of Rxn activity 
showFigure(BRCAProject.comparisons.(compName).functionalComparison.plots.fba.heatmapRxnFluxsum)

% used in terms of metabolite usage 
showFigure(BRCAProject.comparisons.(compName).functionalComparison.plots.fba.heatmapMetsFluxsum)

% how is the cardinality in different subsystems ? 

showFigure(BRCAProject.comparisons.(compName).functionalComparison.plots.fba.heatmapRxnActivityFba)


%% sampling comparison

% how much are metabolites used in a specific subsystem in sampling (average over samples)
showFigure(BRCAProject.comparisons.(compName).samplingComparison.plots.heatmapMetsFluxSum)
% is this trend consistent within the model, or highly variable between samples ? 
showFigure(BRCAProject.comparisons.(compName).samplingComparison.plots.heatmapMetsFluxSumSamples)

% how active are rxns in the different subsystem
showFigure(BRCAProject.comparisons.(compName).samplingComparison.plots.heatmapRxnFluxSum)
% is this trend consistent within the model, or highly variable between samples ? 
showFigure(BRCAProject.comparisons.(compName).samplingComparison.plots.heatmapRxnFluxSumSamples)

%% GENERATING FIGURES ON DEMAND

referenceModel = "consistentMediumConstrainedModel"

% Questions to be answered by plotting data using the functions + getting to know the behaviour of the models in detail
% - How active is a given set of rxns in the fba solution or the sampling solutions ?
% - How much is a given set of metabolites used in the fba solution or the sampling solutions ?
% - How much can the values of a given set of rxns be pushed ? 
% - What is the rxn flux value under the minimal solution for a given set of rxns ? 
% - How do the rxn flux values vary in the sampling? 

% visualize with new functions 

% heatmaps

[rxnsMetId,producingMet,matched] = getRxnIDs(BRCAProject,referenceModel, ["Pentose.*"; "Glycolysis.*"]);

visDiffRxnSetActivitySampling(BRCAProject, compName, rxnsMetId, ["Pentose.* & g6p.*"; "Glycolysis.*"],referenceModel)


visDiffRxnSetActivityFBA(BRCAProject, compName, rxnsMetId, ["Pentose.* & g6p.*"; "Glycolysis.*"],referenceModel)


visDiffMetSetUsageSampling(BRCAProject, compName, rxnsMetId, ["Pentose.* & g6p.*"; "Glycolysis.*"],referenceModel)


visDiffMetSetUsageFBA(BRCAProject, compName, rxnsMetId, ["Pentose.* & g6p.*"; "Glycolysis.*"],referenceModel)


visRxnSetVariability(BRCAProject, compName, rxnsMetId, ["Pentose.* & g6p.*"; "Glycolysis.*"],referenceModel)


visMetSetVariability(BRCAProject, compName, rxnsMetId, ["Pentose.* & g6p.*"; "Glycolysis.*"],referenceModel)

% violins

visSingleMetSamplingDistribution(BRCAProject, compName, rxnsMetId(1), ["Pentose.* & g6p.*"],referenceModel)

visSingleRxnSamplingDistribution(BRCAProject, compName, rxnsMetId(1), ["Pentose.* & g6p.*"],referenceModel,true)

% bar plots


visSingleRxnFBA(BRCAProject, compName, rxnsMetId(1),"FVA",false,"thresholdFlux","none")

%%



visualizeSamplingLandscape(BRCAProject, compName, BRCAProject.models.(referenceModel).model.rxns(rxnsMetId(1)))


