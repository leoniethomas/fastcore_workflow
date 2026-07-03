function checkProjectFormat(project)
% Checks that a project struct conforms to the expected format.
% Only required fields: project.models, at least one model, and model.model.
% Optional fields are validated for type and allowed names if present.

    arguments
        project (1,1) struct
    end

    % Top-level: only 'models' allowed
    actual = fieldnames(project);
    extra = setdiff(actual, {'models'});
    if ~isempty(extra)
        error("Unexpected field(s) in project: %s", strjoin(extra, ", "));
    end

    if ~isfield(project, 'models')
        error("project must have a 'models' field.");
    end

    modelNames = fieldnames(project.models);
    if isempty(modelNames)
        error("project.models must contain at least one model.");
    end

    for i = 1:numel(modelNames)
        name = modelNames{i};
        m = project.models.(name);
        validateModelInProject(m, "project.models." + name);
    end

end

function validateModelInProject(m, path)
% Validate a single model inside project.models

    allowed = {'model', 'sampleMetadata', 'discretizedData', 'expressionData', ...
               'coreReactions', 'settings'};

    actual = fieldnames(m);
    extra = setdiff(actual, allowed);
    if ~isempty(extra)
        error("Unexpected field(s) in %s: %s. Allowed ones are: 'model', 'sampleMetadata', " + ...
              "'discretizedData', 'expressionData', 'coreReactions', 'settings'", ...
              path, strjoin(extra, ", "));
    end

    % Only 'model' is required
    if ~isfield(m, 'model')
        error("Missing required field 'model' in %s.", path);
    end
    mustBeCobraModel(m.model);

    % Type validation for optional fields
    if isfield(m, 'sampleMetadata') && ~istable(m.sampleMetadata)
        error("%s.sampleMetadata must be a table.", path);
    end
    if isfield(m, 'discretizedData') && ~isa(m.discretizedData, 'double')
        error("%s.discretizedData must be a double array.", path);
    end
    if isfield(m, 'expressionData') && ~istable(m.expressionData)
        error("%s.expressionData must be a table.", path);
    end
    if isfield(m, 'coreReactions')
        mustBeVectorOrEmpty(m.coreReactions);
    end

    % Validate settings
    if isfield(m, 'settings')
        validateSettingsInProject(m.settings, path + ".settings", m);
    end

end

function validateSettingsInProject(s, path, m)
% Validate the settings struct

    allowed = {'medium', 'scriptParameters', 'dico', 'objFunction', ...
               'referenceModel', 'mapping', 'optionalSettings'};

    actual = fieldnames(s);
    extra = setdiff(actual, allowed);
    if ~isempty(extra)
        error("Unexpected field(s) in %s: %s. Allowed ones are: 'medium', 'scriptParameters', " + ...
              "'dico', 'objFunction', 'referenceModel', 'mapping', 'optionalSettings'", ...
              path, strjoin(extra, ", "));
    end

    % Type validation
    if isfield(s, 'dico') && ~istable(s.dico)
        error("%s.dico must be a table.", path);
    end
    if isfield(s, 'objFunction') && ~(ischar(s.objFunction) || isstring(s.objFunction))
        error("%s.objFunction must be text.", path);
    end
    if isfield(s, 'referenceModel') && ~(ischar(s.referenceModel) || isstring(s.referenceModel))
        error("%s.referenceModel must be text.", path);
    end
    if isfield(s, 'optionalSettings')
        validateOptionalSettings(s.optionalSettings);
    end

    if isfield(s, 'medium')
        validateMediumInProject(s.medium, path + ".medium");
    end

    if isfield(s, 'scriptParameters')
        validateScriptParametersInProject(s.scriptParameters, path + ".scriptParameters", m);
    end

end

function validateMediumInProject(med, path)
% Validate the medium struct

    allowed = {'mediumComposition', 'manuallySetBoundaries'};

    actual = fieldnames(med);
    extra = setdiff(actual, allowed);
    if ~isempty(extra)
        error("Unexpected field(s) in %s: %s. Allowed ones are: 'mediumComposition', 'manuallySetBoundaries'.", path, strjoin(extra, ", "));
    end

    if isfield(med, 'mediumComposition') && ~istable(med.mediumComposition)
        error("%s.mediumComposition must be a table.", path);
    end
    if isfield(med, 'manuallySetBoundaries')
        validateManuallySetBoundaries(med.manuallySetBoundaries);
    end

end

function validateScriptParametersInProject(sp, path, m)
% Validate scriptParameters and cross-field with sampleMetadata

    allowed = {'sampleLabeling', 'consensusProportion'};

    actual = fieldnames(sp);
    extra = setdiff(actual, allowed);
    if ~isempty(extra)
        error("Unexpected field(s) in %s: %s. Allowed ones are: 'sampleLabeling', 'consensusProportion'.", path, strjoin(extra, ", "));
    end

    % Cross-field validation
    if isfield(sp, 'sampleLabeling')
        if isfield(m, 'sampleMetadata')
            if isempty(m.sampleMetadata) ~= (strlength(sp.sampleLabeling) == 0)
                error("sampleMetadata and sampleLabeling must be provided together or both omitted.");
            end
            if ~isempty(m.sampleMetadata)
                if ~ismember(string(sp.sampleLabeling), string(m.sampleMetadata.Properties.VariableNames))
                    error("'%s' must be a column in sampleMetadata.", sp.sampleLabeling);
                end
            end
        end
    end

    if isfield(sp, 'consensusProportion') && ~isempty(sp.consensusProportion)
        mustBeNonnegative(sp.consensusProportion);
    end

end

%%
% project.models
% project.models.Name1
% project.models.Name1.model
% project.models.Name1.sampleMetadata
% project.models.Name1.discretizedData
% project.models.Name1.expressionData
% project.models.Name1.coreReactions
% project.models.Name1.settings
% project.models.Name1.settings.medium
% project.models.Name1.settings.medium.mediumComposition
% project.models.Name1.settings.medium.manuallySetBoundaries
% project.models.Name1.settings.medium.manuallySetBoundaries.closedImports
% project.models.Name1.settings.medium.manuallySetBoundaries.closedExports
% project.models.Name1.settings.medium.manuallySetBoundaries.unconstrainedImports
% project.models.Name1.settings.medium.manuallySetBoundaries.unconstrainedExports
% project.models.Name1.settings.scriptParameters
% project.models.Name1.settings.scriptParameters.consensusProportion
% project.models.Name1.settings.scriptParameters.sampleLabeling
% project.models.Name1.settings.dico
% project.models.Name1.settings.objFunction
% project.models.Name1.settings.referenceModel
% project.models.Name1.settings.mapping
% project.models.Name1.settings.optionalSettings
% project.models.Name1.settings.optionalSettings.medium
% project.models.Name1.settings.optionalSettings.notMediumConstrained
% project.models.Name1.settings.optionalSettings.func

% project.models.Name1.analysis
% project.models.Name1.analysis.analysis_id
% project.models.Name1.analysis.analysis_id.parameters
% project.models.Name1.analysis.analysis_id.FBA
% project.models.Name1.analysis.analysis_id.FVA
% project.models.Name1.analysis.analysis_id.loopStatus
% project.models.Name1.analysis.analysis_id.singleGeneDeletion
% project.models.Name1.analysis.analysis_id.doubleGeneDeletion
% project.models.Name1.analysis.analysis_id.sampling
% project.models.Name1.analysis.analysis_id.loopless
% project.models.Name1.analysis.analysis_id.kld