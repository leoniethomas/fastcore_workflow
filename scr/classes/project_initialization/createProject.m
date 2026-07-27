function project = createProject(params)
% this function creates a project object ready to be used in the pipeline
% for analysis

% all fields 
% - modelName
% - contextSpecificModel (rFastcormics)
% - expressionData
% - discretized (rFastcormics)
% - geneNames (rFastcormics)
% - dico (rFastcormics)
% - objFunction (rFastcormics)
% - consensusProportion (rFastcormics)
% - optionalSettings: (rFastcormics)
%      - medium
%      - notMediumConstrained
%      - .func
% - referenceModel (rFastcormics)
% - mapping (rFastcormics)
% - coreReactions (rFastcormics)
% - mediumComposition
% - manuallySetBoundaries:
%      - closedImport
%      - closedExport
%      - openedImport
%      - openedExport
% - sampleMetadata
% - sampleLabeling

arguments
    params (1,:) cell
end

% Initialize project
project = struct();
project.models = struct();

% Loop through models
for i = 1:numel(params)
    paramsForModel = params{i};
    
    % Validate each struct's fields
    paramsForModel = validateParamsForPipeline(paramsForModel);
    
    % Initiate a struct per model
    project.models.(paramsForModel.modelName) = formatParamsForModel(paramsForModel);
end

end
