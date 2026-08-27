function [fluxSet,figs] = visSingleRxnSamplingDistribution(project, compName, rxnSet,rxnSetLabel,referenceModel,addKLDValues)
arguments
    project 
    compName (1,1) string
    rxnSet (1,:) cell
    rxnSetLabel (1,:) string
    referenceModel (1,1) string
    addKLDValues (1,1) logical = false
end

% 
% This function generates a violin that shows the distribution of flux
% values overall samples for one rxns. One violin plot per defined index in
% rxnSet is generated.
% USAGE: 
%       visSingleRxnSamplingDistribution(project, compName, rxnSet,rxnSetLabel,referenceModel) 
% INPUT: 
%       project:        object generated through the pipeline
%                       (singleModelAnalysis and modelsComparison have to be run before)
%       compName:       the comparison to visualize in this figure. The names
%                       of the comparisons already present in the object can be seen using
%                       project.comparisons
%       rxnSet:        The sets of rxns to be visualized 
%       rxnSetLabel:   The Labels that should be used in the figure
%       referenceModel: The name of the reference model used.
%   optional: 
%       addKLDValue:    Specify whether you want the significance of
%                       dissimilarity between the sampling distribution to
%                       be displayed in the violin plot. (default:false)
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

% check that the reference model used when retrieving the rxnIDs is the
% same as in the comparison, otherwise the ids are wrongly assigned

if project.comparisons.(compName).referenceModel ~= referenceModel
    error('The reference Model you gave as input (%s) is not the same reference Model defined in the specified comparison object (%s). Check the spelling! The reference model when retrieving the rxnIDS with getRxnIDs needs to be the same as the one specified here! Otherwise the idx are misaligned.  "%s" keeps referenceModel = "%s".Choose a different identifier to create a new comparison.', ...
          referenceModel, project.comparisons.(compName).referenceModel);
end

%

%
rxnLabel    = matlab.lang.makeValidName(rxnSetLabel);  % used as field name

[fluxSet,figs] = visualizeFlux(project, compName,rxnSet,rxnLabel,"orderedSamples","all",'on',addKLDValues)

end