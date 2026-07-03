function validModels = checkStructureForSingleModelAnalysis(project, modelList)
% Validates that proj is a structure containing a field "models"
% which is a structure containing model objects listed in modelList.
% If modelList is empty, checks all fields in proj.models.
% Throws an error and stops if a problem is detected.
% Returns a cell array of valid model names.

    % Check proj is a structure
    if ~isstruct(project)
        error('Invalid input: proj must be a structure.');
    end

    % Check "models" field exists
    if ~isfield(project, 'models')
        error(['Invalid structure: missing field "models". ' ...
               'Expected proj.models to be a structure containing structures named as in modelList.']);
    end

    % Check models is a structure
    if ~isstruct(project.models)
        error('Invalid structure: field "models" must be a structure.');
    end

    % If modelList is empty, use all fields in proj.models
    if isempty(modelList) || (ischar(modelList) && isempty(modelList))
        modelList = fieldnames(project.models);
        if isempty(modelList)
            error('No models found in proj.models.');
        end
    end

    % Check each model in modelList
    missingModels = {};
    invalidModels = {};
    missingModelField = {};
    validModels = {};

    for i = 1:numel(modelList)
        name = modelList{i};

        % Check model struct exists
        if ~isfield(project.models, name)
            missingModels{end+1} = name;
            continue
        end

        modelStruct = project.models.(name);

        % Check it is a structure
        if ~isstruct(modelStruct)
            invalidModels{end+1} = name;
            continue
        end

        % Check "model" field exists
        if ~isfield(modelStruct, 'model')
            missingModelField{end+1} = name;
            continue
        end

        % Check type of model
        modelObj = modelStruct.model;
        if ~isstruct(modelObj) || ...
           ~all(isfield(modelObj, {'S','rxns','mets','lb','ub', 'genes', 'subSystems'}))
            invalidModels{end+1} = name;
            continue
        end

        % Model is valid
        validModels{end+1} = name;
    end

    % Report results
    if ~isempty(missingModels) || ~isempty(invalidModels) || ~isempty(missingModelField)

        if ~isempty(missingModels)
            fprintf('Missing model structures in proj.models:');
            fprintf(' - %s', missingModels{:});
        end

        if ~isempty(missingModelField)
            fprintf('Missing "model" field in model structures:');
            fprintf(' - %s', missingModelField{:});
        end

        if ~isempty(invalidModels)
            fprintf('Invalid model objects (not class "model"):');
            fprintf(' - %s', invalidModels{:});
        end

        error('Structure validation failed.');
    end

    disp('Structure OK: all models are present and valid.');

end
