function [structureOk] = checkStructureForSingleModelAnalysis(proj, modelList)
% Validates that proj is a structure containing a field "models"
% which is a structure containing model objects listed in modelList

% add something for checking by default all the models in models ?

    structureOk = false;

    % Check proj is a structure
    if ~isstruct(proj)
        error('Invalid input: proj must be a structure.');
        return
    end

    % Check "models" field exists
    if ~isfield(proj, 'models')
        error(['Invalid structure: missing field "models". ' ...
               'Expected proj.models to be a structure containing structures named as in modelList.']);
        return
    end

    % Check models is a structure
    if ~isstruct(proj.models)
        error(['Invalid structure: field "models" must be a structure.']);
        return
    end

    % Check each model in modelList
    missingModels = {};
    invalidModels = {};
    missingModelField = {};

    for i = 1:numel(modelList)
        name = modelList{i};

        % Check model struct exists
        if ~isfield(proj.models, name)
            missingModels{end+1} = name;
            continue
        end

        modelStruct = proj.models.(name);

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
           ~all(isfield(modelObj, {'S','rxns','mets','lb','ub', 'genes', 'subSystems'})) % see if other fields are required
            invalidModels{end+1} = name;
        end
    end

    % Report results
    if ~isempty(missingModels) || ~isempty(invalidModels) || ~isempty(missingModelField)

        if ~isempty(missingModels)
            fprintf('Missing model structures in proj.models:\n');
            fprintf('  - %s\n', missingModels{:});
        end

        if ~isempty(missingModelField)
            fprintf('Missing "model" field in model structures:\n');
            fprintf('  - %s\n', missingModelField{:});
        end

        if ~isempty(invalidModels)
            fprintf('Invalid model objects (not class "model"):\n');
            fprintf('  - %s\n', invalidModels{:});
        end

        structureOk = false;
        return
    end

    % Success
    disp('Structure OK: all models are present and valid.');
    structureOk = true;

end

