function project = addModelsToProject(project, params)
% this function adds one or several model to an already existing project.

% all fields 
% - modelName
% - contextSpecificModel (rFastcormics)
% - expressionData
% - discretized (rFastcormics)
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
    project
    params (1,:) cell
end

% Check whether the project is in the correct format
fprintf("Checking project format before adding model. \n")
checkProjectFormat(project);

% Add extra models one by one
for i = 1:numel(params)
    paramsForModel = params{i};
    
    % Validate each struct's fields
    nv = struct2nv(paramsForModel); % converts in name-value pairs
    paramsForModel = validateParamsForPipeline(nv{:});

    % Initiate a struct per model
    if isfield(project.models, paramsForModel.modelName)
        answer = input("Model '" + paramsForModel.modelName + "' already exists. Overwrite? [y/n]: ", 's');
        if ~strcmpi(answer, 'y')
            fprintf("Model '" + paramsForModel.modelName + "' skipped.")
            continue
        end
    end

    project.models.(paramsForModel.modelName) = formatParamsForModel(paramsForModel);

end

end

function args = struct2nv(s)
    fn = fieldnames(s);
    vals = struct2cell(s);
    args = reshape([fn, vals]', 1, []);
end
