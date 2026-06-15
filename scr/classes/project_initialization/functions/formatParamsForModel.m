function modelShortcut = formatParamsForModel(paramsForModel)
%   Reformats pipeline parameters into a model struct.
%   Takes a paramsForModel structure (as validated by validateParamsForPipeline)
%   and reorganizes its fields into the nested structure format expected by the
%   pipeline, separating model data, medium settings, and script parameters.

modelShortcut = struct();
modelShortcut.model = paramsForModel.contextSpecificModel;

% Optional fields
if length(fieldnames(paramsForModel)) > 2
    modelShortcut.settings = struct();
    modelShortcut.settings.medium = struct();
    modelShortcut.settings.scriptParameters = struct();
    modelShortcut.settings.medium.manuallySetBoundaries = struct();
    
    % Optional fields
    if isfield(paramsForModel, 'sampleMetadata')
        modelShortcut.sampleMetadata = paramsForModel.sampleMetadata;
    end

    if isfield(paramsForModel, 'discretizedData')
        modelShortcut.discretizedData = paramsForModel.discretizedData;
    end

    if isfield(paramsForModel, 'expressionData')
        modelShortcut.expressionData = paramsForModel.expressionData;
    end
    
    if isfield(paramsForModel, 'coreReactions')
        modelShortcut.coreReactions = paramsForModel.coreReactions;
    end
    
    if isfield(paramsForModel, 'dico')
        modelShortcut.settings.dico = paramsForModel.dico;
    end

    if isfield(paramsForModel, 'mediumComposition')
        modelShortcut.settings.medium.mediumComposition = paramsForModel.mediumComposition;
    end
    
    if isfield(paramsForModel, 'manuallySetBoundaries')
        modelShortcut.settings.medium.manuallySetBoundaries = paramsForModel.manuallySetBoundaries;
        if isfield(paramsForModel, 'optionalSettings')
            if isfield(paramsForModel.optionalSettings, 'notMediumConstrained')
                modelShortcut.settings.medium.manuallySetBoundaries.unconstrainedImports = unique([paramsForModel.manuallySetBoundaries.unconstrainedImports, ...
                    paramsForModel.optionalSettings.notMediumConstrained]);
            end
        end
    end
    
    if isfield(paramsForModel, 'consensusProportion')
        modelShortcut.settings.scriptParameters.consensusProportion = paramsForModel.consensusProportion;
    end
    
    if isfield(paramsForModel, 'optionalSettings')
        modelShortcut.settings.optionalSettings = paramsForModel.optionalSettings;
        if isfield(paramsForModel.optionalSettings, 'notMediumConstrained') && ~isfield(paramsForModel.manuallySetBoundaries, 'unconstrainedImports')
            modelShortcut.settings.medium.manuallySetBoundaries.unconstrainedImports = paramsForModel.optionalSettings.notMediumConstrained;
        end
    end
    
    if isfield(paramsForModel, 'objFunction')
        modelShortcut.settings.objFunction = paramsForModel.objFunction;
    end
    
    if isfield(paramsForModel, 'referenceModel')
        modelShortcut.settings.referenceModel = paramsForModel.referenceModel;
    end
    
    if isfield(paramsForModel, 'mapping')
        modelShortcut.settings.mapping = paramsForModel.mapping;
    end

    if isfield(paramsForModel, 'sampleLabeling')
        modelShortcut.settings.scriptParameters.sampleLabeling = paramsForModel.sampleLabeling;
    end

    % Remove empty fields (except in .model fields)
    modelShortcut = removeRecursivelyEmptyFields(modelShortcut, {'model'});

end

end