function [fluxSet,figs] = visSingleMetSamplingDistribution(project, compName, rxnSet,rxnSetLabel,referenceModel)
arguments
    project 
    compName (1,1) string
    rxnSet (1,1) cell
    rxnSetLabel (1,1) string
    referenceModel (1,1) string
end

% 
% This function generates one violin plot that shows the fluxsum value distribution for
% all the metabolites taking part in the defined rxns. Each dot in the
% violin plot is one sample.
% USAGE: 
%       visSingleMetSamplingDistribution(project, compName, rxnSet,rxnSetLabel,referenceModel) 
% INPUT: 
%       project:        object generated through the pipeline
%                       (singleModelAnalysis and modelsComparison have to be run before)
%       compName:       the comparison to visualize in this figure. The names
%                       of the comparisons already present in the object can be seen using
%                       project.comparisons
%       rxnSet:        One set of rxns to be visualized 
%       rxnSetLabel:   The label to be displayed for this rxnSet in the plot.
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


% check that the reference model used when retrieving the rxnIDs is the
% same as in the comparison, otherwise the ids are wrongly assigned

if project.comparisons.(compName).referenceModel ~= referenceModel
    error('The reference Model you gave as input (%s) is not the same reference Model defined in the specified comparison object (%s). Check the spelling! The reference model when retrieving the rxnIDS with getRxnIDs needs to be the same as the one specified here! Otherwise the idx are misaligned.  "%s" keeps referenceModel = "%s".Choose a different identifier to create a new comparison.', ...
          referenceModel, project.comparisons.(compName).referenceModel);
end

%
rxnLabel    = matlab.lang.makeValidName(rxnSetLabel);  % used as field name

[fluxSet,figs] = visualizeFluxsum(project, compName,[],rxnSet,rxnLabel,"violin",...
                                       false,false,"orderedSamples","outgoing",referenceModel,'on')

end