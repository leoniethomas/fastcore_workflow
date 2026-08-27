function [fluxSet,figs] = visDiffRxnSetActivitySampling(project, compName, rxnSets,rxnSetLabels,referenceModel,metIdx)
arguments
    project 
    compName (1,1) string
    rxnSets (1,:) cell
    rxnSetLabels (1,:) string
    referenceModel (1,1) string
end

% 
% This function generates a heatmap that shows how active the defined rxns
% are, by building the sum of all the values for those reactions in one
% sample. The resulting heatmap entail one value for a RxnSet for one sample.
% This is calculated based on all solutions from the performed sampling,
% and an average is build over all samples to get one value per rxnSet in
% the heatmap.
% The usage here is defined as the flux sum of all the rxns that consume
% this metabolite in the defined rxnSet. 
% USAGE: 
%       visDiffMetSetUsageSampling(project, compName, rxnSets,rxnSetLabels,referenceModel) 
% INPUT: 
%       project:        object generated through the pipeline
%                       (singleModelAnalysis and modelsComparison have to be run before)
%       compName:       the comparison to visualize in this figure. The names
%                       of the comparisons already present in the object can be seen using
%                       project.comparisons
%       rxnSets:        The sets of rxns to be visualized 
%       rxnSetLabels:   The Labels that should be used in the figure
%       referenceModel: The name of the reference model used.
%
% OUTPUT:
%       figs:           Figure object generated and displayed by the
%                       function
%
%
% Code Example: 
% after running singleModelAnalysis and modelsComparison this code can be
% used: 
% 
% [rxnsMetId,producingMet,matched] = getRxnIDs(BRCAProject,referenceModel, ["Pentose.* & g6p.*"; "Glycolysis.*"]);
% visDiffMetSetUsageFBA(BRCAProject, compName, rxnsMetId, ["Pentose.* & g6p.*"; "Glycolysis.*"],referenceModel)
%

% check project format: 
% comparison specified is there ? has the required fields ? 

if ~isfield(project, 'comparisons')
    error('Project object does not contain a comparison object. After running the singleModels analysis you still have to run modelsComparison before beeing able to use this function!')
end
if ~isfield(project.comparisons, compName)
    error('The comparison name you gave as an input is not available in the object. Check your spelling!')
end
if ~isfield(project.comparisons.(compName),'samplingComparison')
    error('The comparison object you specified does not entail a samplingComparison. Run modelsComparison(project, modelsToCompare, referenceModel, compID, ["samplingComparison"]).')
end
if ~isfield(project.comparisons.(compName).samplingComparison, 'orderedSamples')
    error('The comparisons object does not entail a valid samplingComparison object rerun the functionalComparison: modelsComparison(project, modelsToCompare, referenceModel, compID, ["samplingComparison"])')
end

% rxnSets and rxSetLabels need to be the same length
if length(rxnSets) ~= length(rxnSetLabels)
    error('The rxnSets and the rxnSetLabels you specified do not have the same length, you need to give as many Lables as rxnID sets. ')
end

% check that the reference model used when retrieving the rxnIDs is the
% same as in the comparison, otherwise the ids are wrongly assigned

if project.comparisons.(compName).referenceModel ~= referenceModel
    error('The reference Model you gave as input (%s) is not the same reference Model defined in the specified comparison object (%s). Check the spelling! The reference model when retrieving the rxnIDS with getRxnIDs needs to be the same as the one specified here! Otherwise the idx are misaligned.  "%s" keeps referenceModel = "%s".Choose a different identifier to create a new comparison.', ...
          referenceModel, project.comparisons.(compName).referenceModel);
end

%

[fluxSet,figs] = visualizeFluxsum(project, compName,metIdx,rxnSets,rxnSetLabels,"heatmap",...
                                       false,false,"orderedllSamples","reactions",referenceModel,'on')

end