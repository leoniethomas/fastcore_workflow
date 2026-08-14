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
load('BRCAProjectNewFormat_0508.mat');
% load('defaultParametersAnalysis.csv'); % can be modified according to the
% user % this does not work for me, throws an error 
parametersAnalysis = readtable('defaultParametersAnalysis.csv', 'TextType', 'string');

%% COMPUTING THE COMPARISON
modelsToCompare = {"Control", "StageIV", "StageI"};
[BRCAProject, analysisIDs] = chooseActiveAnalysis(BRCAProjectBis, modelsToCompare); % by default, the most recent analyses will be chosen 

% analysisToChoose = {'analysis_20260717_1704', 'analysis_20260717_1757', 'analysis_20260717_1841'};
% [BRCAProject, analysisIDs] = chooseActiveAnalysis(BRCAProject, modelsToCompare, analysisToChoose);

%% VISUALIZING AUTOMATICALLY GENERATED FIGURES

referenceModel = "consistentMediumConstrainedModel"; 
comparisonList = ["structuralComparison",...
    "functionalComparison", ...
    "samplingComparison"];
compID = "tutorialComparison_20260728aa"; 

%[BRCAProjectBis, comparisonName] = modelsComparison(BRCAProjectBis, modelsToCompare, referenceModel, compID); % only structural comp
[BRCAProject, comparisonName] = modelsComparison(BRCAProject, modelsToCompare, referenceModel, compID, comparisonList);

%% Showing the Default figures generated during comparison 


% --- Structural Analysis Plots

% how similar are the models to one another 
% jaccard similarity over the model structure, meaning are rxns included in
% the model
showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.structuralComparison.plots.jaccardDist.rxns)
% most reactions are shared between the models, which is to be expected 
% what can also be observed is that the cancer models are more similar than
% to the normal model

% lets check out the size and overlap of the models, in terms of genes,
% metabolites and reactions
showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.structuralComparison.plots.intersections.genes)
showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.structuralComparison.plots.intersections.mets)
showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.structuralComparison.plots.intersections.rxns)
% in these venn diagramms we see in absolute numbers what is indicated by
% the jaccard similarity heatmap

% how were the genes in the model discretized ? 
showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.structuralComparison.plots.dataDiscretization)
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
showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.structuralComparison.plots.coreReactions)
% the next question then is how many of those core reactions are specific
% or are the all shared ? 
% The models share a large amount of core reactions, which makes sense but we also see that the control samples have the highest number of specific rxns, 
% while the overlab between the two cancer models are larger than between
% the cancer models and control, stage I has the least number of specific
% reactions

% where do those core reactions come from, pathway wise 
showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.structuralComparison.plots.coreReactionsIntersections)
% the important observation here is that the differences in core reaction
% is not only due to Transport or exchangers, but there is also other rxns
% that make up the model specific core reactions



%the next question is where do the differences between the models come from
%in terms of rxns presecence, so structural difference
showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.structuralComparison.plots.reactionPathwayPresence)


%%
showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.functionalComparison.plots.import)

showFigure(BRCAProject.comparisons.Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa.functionalComparison.plots.export)


%% GENERATING FIGURES ON DEMAND

referenceModel = "consistentMediumConstrainedModel"
[rxnsMetId,producingMet,matched] = getRxnIDs(BRCAProject,referenceModel, "Glycolysis.* & g6p[c");
[fluxsumSets,plot1] = visualizeFlux(BRCAProject, "Control_vs_StageI_vs_StageIV__tutorialComparison_20260728aa",{rxnsMetId },...
                                                                     ["Lactate and glucos glutamine export"],"all");


%%

