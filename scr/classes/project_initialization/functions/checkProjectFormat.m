function checkProjectFormat(project, modelList, analysisList)
% Checks that a project struct conforms to the expected format.
% Only required fields: project.models, at least one model, and model.model.
% Optional fields are validated for type and allowed names if present.
%
% USAGE:
%   checkProjectFormat(project)
%   checkProjectFormat(project, modelList)
%   checkProjectFormat(project, modelList, analysisList)
%
% If modelList is empty, all models are validated.
% If modelList is provided, only those models are validated (must exist).
% If analysisList is provided, modelList must also be provided and have
% the same length. Only the specified analyses are validated, each
% analysis must exist in its corresponding model (same order).

    arguments
        project (1,1) struct
        modelList (1,:) string = strings(1,0)
        analysisList (1,:) string = strings(1,0)
    end

    % Top-level: only 'models' and 'comparisons' allowed
    actual = fieldnames(project);
    extra = setdiff(actual, {'models', 'comparisons'});
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

    % Validate analysisList constraints
    if ~isempty(analysisList)
        if isempty(modelList)
            error("analysisList cannot be provided without modelList.");
        end
        if numel(analysisList) ~= numel(modelList)
            error("analysisList must have the same number of elements as modelList " + ...
                "(got %d analyses for %d models).", numel(analysisList), numel(modelList));
        end
    end

    % Determine which models to validate
    if isempty(modelList)
        modelsToValidate = string(modelNames);
    else
        modelsToValidate = modelList;
        availableModels = string(modelNames);
        missing = setdiff(modelList, availableModels);
        if ~isempty(missing)
            error("Requested model(s) not found in project.models: %s", strjoin(missing, ", "));
        end
    end

    % Loop over models
    for i = 1:numel(modelsToValidate)
        name = modelsToValidate(i);
        m = project.models.(name);

        if ~isempty(analysisList)
            analysisName = analysisList(i);
            validateModelInProject(m, "project.models." + name, analysisName);
        else
            validateModelInProject(m, "project.models." + name);
        end
    end

end

function validateModelInProject(m, path, analysisName)
% Validate a single model inside project.models.
% If analysisName is provided and non-empty, only that analysis entry
% is validated. Otherwise, all analyses are validated.

    allowed = {'model', 'sampleMetadata', 'discretizedData', 'expressionData', ...
        'coreReactions', 'mappedDiscretizedRxnsAllSamples', 'mappedDiscretizedRxns', ...
        'settings', 'analysis'};

    actual = fieldnames(m);
    extra = setdiff(actual, allowed);
    if ~isempty(extra)
        error("Unexpected field(s) in %s: %s. Allowed ones are: 'model', 'sampleMetadata', " + ...
            "'discretizedData', 'expressionData', 'coreReactions', 'mappedDiscretizedRxnsAllSamples', 'mappedDiscretizedRxns', 'settings', 'analysis'", ...
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
    
    if isfield(m, 'discretizedData')
        if ~istable(m.discretizedData)
            error("%s.discretizedData must be a table.", path);
        end
        requiredCols = {'geneIds', 'value'};
        missingCols = setdiff(requiredCols, string(m.discretizedData.Properties.VariableNames));
        if ~isempty(missingCols)
            error("%s.discretizedData must contain columns: %s. Missing: %s.", ...
                path, strjoin(requiredCols, ", "), strjoin(missingCols, ", "));
        end

        % Checking of m.mappingDiscretizedRxnsAllSamples (int8) and
        % m.mappedDiscretizedRxns (int8)
        % Same nb of rows as discretizedData
        nRows = size(m.model.rxns, 1);
    
        % --- mappingDiscretizedRxnsAllSamples : int8 ---
        if ~isa(m.mappedDiscretizedRxnsAllSamples, 'int8')
            error("%s.mappedDiscretizedRxnsAllSamples must be int8 (current: %s).", ...
                  path, class(m.mappedDiscretizedRxnsAllSamples));
        end
        if size(m.mappedDiscretizedRxnsAllSamples, 1) ~= nRows
            name = string(regexp(path, '[^.]+$', 'match', 'once'));
            error("%s.mappedDiscretizedRxnsAllSamples must have %d rows (same as number of reactions in %s model), has %d.", ...
                  path, nRows, name, size(m.model.rxns, 1));
        end
    
        % --- mappedDiscretizedRxns : int8 ---
        if ~isa(m.mappedDiscretizedRxns, 'int8')
            error("%s.mappedDiscretizedRxns must be int8 (current: %s).", ...
                  path, class(m.mappedDiscretizedRxns));
        end
        if size(m.mappedDiscretizedRxns, 1) ~= nRows
            name = string(regexp(path, '[^.]+$', 'match', 'once'));
            error("%s.mappedDiscretizedRxns must have %d rows (same as number of reactions in %s model), has %d.", ...
                  path, nRows, name, size(m.model.rxns, 1));
        end
    end
    
    if isfield(m, 'expressionData')
        if ~istable(m.expressionData)
            error("%s.expressionData must be a table.", path);
        end
        requiredCols = {'geneIds', 'expression'};
        missingCols = setdiff(requiredCols, string(m.expressionData.Properties.VariableNames));
        if ~isempty(missingCols)
            error("%s.expressionData must contain columns: type: %s. Missing: %s.", ...
                path, strjoin(requiredCols, ", "), strjoin(missingCols, ", "));
        end
    end
    
    if isfield(m, 'coreReactions')
        mustBeVectorOrEmpty(m.coreReactions);
    end

    % Validate settings
    if isfield(m, 'settings')
        validateSettingsInProject(m.settings, path + ".settings", m);
    end

     % === Cross-field validation: dico, discretizedData, expressionData, model.genes ===
     hasDiscretized = isfield(m, 'discretizedData');
     hasExpression = isfield(m, 'expressionData');
     hasSettings = isfield(m, 'settings');
     hasDico = hasSettings && isfield(m.settings, 'dico');
     nbModelGenes = numel(m.model.genes);
     modelGenes = regexprep(string(m.model.genes), "\.[0-9]+$", "");
    
     % If discretizedData or expressionData -> dico required with geneIdsInData column
     if (hasDiscretized || hasExpression) && ~hasDico
        error("%s: dico is required in settings when discretizedData or expressionData is present.", path);
     end
    
     if hasDico
     % dico must have geneIdsInData if discretizedData or expressionData
        if (hasDiscretized || hasExpression)
            if ~ismember('geneIdsInData', string(m.settings.dico.Properties.VariableNames))
                error("%s.settings.dico must contain column 'geneIdsInData' when discretizedData or expressionData is present.", path);
            end
        end
    
        % dico must have same number of rows as model.genes
        if height(m.settings.dico) ~= nbModelGenes
            error("%s.settings.dico must have %d rows (same as model.genes), got %d.", ...
                path, nbModelGenes, height(m.settings.dico));
        end
    
        % dico.geneIdsInModel must be in same order as model.genes
        dicoModelGenes = regexprep(string(m.settings.dico.geneIdsInModel), "\.[0-9]+$", "");
        validMask = ~isempty(dicoModelGenes);
        if ~isequal(dicoModelGenes(validMask), modelGenes(validMask))
            error("%s.settings.dico.geneIdsInModel must be in the same order as model.genes.", path);
        end
     end
    
     % discretizedData: same number of rows as dico and model.genes, same gene order as dico.geneIdsInData
     if hasDiscretized
        if height(m.discretizedData) ~= nbModelGenes
            error("%s.discretizedData must have %d rows (same as model.genes and dico), got %d.", ...
                path, nbModelGenes, height(m.discretizedData));
        end
    
        % geneIds in discretizedData must be in same order as dico.geneIdsInData
        discGenes = regexprep(string(m.discretizedData.geneIds), "\.[0-9]+$", "");
        dicoDataGenes = regexprep(string(m.settings.dico.geneIdsInData), "\.[0-9]+$", "");
         if ~isequal(discGenes, dicoDataGenes)
            error("%s.discretizedData.geneIds must be in the same order as dico.geneIdsInData.", path);
         end
     end

    % Validate analysis
    if nargin >= 3 && strlength(analysisName) > 0
        if ~isfield(m, 'analysis')
            error("Model %s does not have an 'analysis' field, but analysis '%s' was requested.", ...
                path, analysisName);
        end
        validateAnalysisInProject(m.analysis, path + ".analysis", m.model, analysisName);
    else
        if isfield(m, 'analysis')
            validateAnalysisInProject(m.analysis, path + ".analysis", m.model);
        end
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

    if isfield(s, 'dico')
        if ~istable(s.dico)
            error("%s.dico must be a table.", path);
        end
        requiredDicoCols = {'geneIdsInModel'};
        missingCols = setdiff(requiredDicoCols, string(s.dico.Properties.VariableNames));
        if ~isempty(missingCols)
            error("%s.dico must contain at least column: %s. Missing: %s.", ...
                path, strjoin(requiredDicoCols, ", "), strjoin(missingCols, ", "));
        end
    end

    if isfield(s, 'objFunction') && ~(ischar(s.objFunction) || isstring(s.objFunction))
        error("%s.objFunction must be text.", path);
    end

    if isfield(s, 'referenceModel') && ~(ischar(s.referenceModel) || isstring(s.referenceModel))
        error("%s.referenceModel must be text.", path);
    end

    if isfield(s, 'mapping')
        if ~(issparse(s.mapping) && isa(s.mapping, 'double'))
            error("s.mapping must be a sparse double matrix.");
        end
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

function validateAnalysisInProject(a, path, model, analysisName)
% Validate the analysis struct of a model.
% All fields in analysis (including 'active') are analysis entries.
% The 'active' entry has an additional 'analysisId' subfield in each
% struct subfield (FBA, FVA, singleGeneDeletion, etc.).
%
% If analysisName is provided and non-empty, only that entry is validated.
% Otherwise, all entries are validated.

    if ~isstruct(a)
        error("%s must be a struct.", path);
    end

    nbRxns = numel(model.rxns);

    if nargin >= 4 && strlength(analysisName) > 0
        if ~isfield(a, char(analysisName))
            error("Analysis '%s' not found in %s.", analysisName, path);
        end
        validateAnalysisEntry(a.(char(analysisName)), path + "." + analysisName, nbRxns, char(analysisName));
    else
        entryNames = fieldnames(a);
        for i = 1:numel(entryNames)
            name = entryNames{i};
            validateAnalysisEntry(a.(name), path + "." + name, nbRxns, name);
        end
    end

end

function validateAnalysisEntry(entry, path, nbRxns, entryName)
% Validate a single analysis entry (active or analysis_id).
% All sub-fields are optional, but if present they must match the expected format.
% The 'active' entry has an additional 'analysisId' (char) subfield in each
% struct subfield (FBA, FVA, singleGeneDeletion, doubleGeneDeletion, sampling, kld).

    isActive = strcmp(entryName, 'active');

    allowed = {'parameters', 'FBA', 'FVA', 'singleGeneDeletion', ...
        'doubleGeneDeletion', 'sampling', 'kld'};

    actual = fieldnames(entry);
    extra = setdiff(actual, allowed);
    if ~isempty(extra)
        error("Unexpected field(s) in %s: %s. Allowed ones are: 'parameters', 'FBA', 'FVA', " + ...
            "'singleGeneDeletion', 'doubleGeneDeletion', 'sampling', 'kld'.", ...
            path, strjoin(extra, ", "));
    end

    % Determine number of samples from parameters table
    nSamples = [];
    if isfield(entry, 'parameters')
        if ~istable(entry.parameters)
            warning("%s.parameters should be a table with following columns: 'Parameter', 'Analysis', 'Value'.", path);
        end
        requiredCols = {'Parameter', 'Analysis', 'Value'};
        missingCols = setdiff(requiredCols, entry.parameters.Properties.VariableNames);
        if ~isempty(missingCols)
            warning("%s.parameters should be a table with following columns: 'Parameter', 'Analysis', 'Value'.", path);
        end
        nSamples = getNumberOfSamples(entry.parameters);
    end

    % Validate FBA (struct; analysisId subfield only in active)
    if isfield(entry, 'FBA')
        if ~isstruct(entry.FBA)
            error("%s.FBA must be a struct.", path);
        end
        if isActive
            allowedFba = {'analysisId', 'basis', 'f', 'f0', 'f1', 'f2', 'lpmethod', 'obj', 'origStat', 'origStatText', 'solver', 'stat', 'time', 'v', 'vars_v', 'x'};
            actualFba = fieldnames(entry.FBA);
            extraFba = setdiff(actualFba, allowedFba);
            if ~isempty(extraFba)
                error("Unexpected field(s) in %s.FBA: %s. Allowed ones are: 'analysisId'.", ...
                    path, strjoin(extraFba, ", "));
            end
            if isfield(entry.FBA, 'analysisId') && ~(ischar(entry.FBA.analysisId) || isstring(entry.FBA.analysisId))
                error("%s.FBA.analysisId must be a char or string.", path);
            end
        end
    end

    % Validate FVA (struct with minMaxFluxes, loopStatus, and analysisId if active)
    if isfield(entry, 'FVA')
        if ~isstruct(entry.FVA)
            error("%s.FVA must be a struct.", path);
        end
        if isActive
            allowedFva = {'minMaxFluxes', 'loopStatus', 'analysisId'};
        else
            allowedFva = {'minMaxFluxes', 'loopStatus'};
        end
        actualFva = fieldnames(entry.FVA);
        extraFva = setdiff(actualFva, allowedFva);
        if ~isempty(extraFva)
            error("Unexpected field(s) in %s.FVA: %s. Allowed ones are: %s.", ...
                path, strjoin(extraFva, ", "), strjoin(allowedFva, ", "));
        end

        % Validate minMaxFluxes (table with minFlux, maxFlux, rows = nb of rxns)
        if isfield(entry.FVA, 'minMaxFluxes')
            if ~istable(entry.FVA.minMaxFluxes)
                error("%s.FVA.minMaxFluxes must be a table.", path);
            end
            if ~all(ismember({'minFlux', 'maxFlux'}, entry.FVA.minMaxFluxes.Properties.VariableNames))
                error("%s.FVA.minMaxFluxes must have 'minFlux' and 'maxFlux' columns.", path);
            end
            if size(entry.FVA.minMaxFluxes, 1) ~= nbRxns
                error("%s.FVA.minMaxFluxes must have %d rows (one per reaction).", path, nbRxns);
            end
        end

        % Validate loopStatus (logical, nbRxns x 1)
        if isfield(entry.FVA, 'loopStatus')
            if ~islogical(entry.FVA.loopStatus)
                error("%s.FVA.loopStatus must be logical.", path);
            end
            if ~isequal(size(entry.FVA.loopStatus), [nbRxns, 1])
                error("%s.FVA.loopStatus must be a %dx1 logical vector.", path, nbRxns);
            end
        end

        % Validate analysisId (char, only in active)
        if isActive && isfield(entry.FVA, 'analysisId')
            if ~(ischar(entry.FVA.analysisId) || isstring(entry.FVA.analysisId))
                error("%s.FVA.analysisId must be a char or string.", path);
            end
        end
    end

    % Validate singleGeneDeletion (struct with sub-fields + analysisId if active)
    if isfield(entry, 'singleGeneDeletion')
        if ~isstruct(entry.singleGeneDeletion)
            error("%s.singleGeneDeletion must be a struct.", path);
        end
        if isActive
            allowedSgd = {'analysisId', 'grRatio', 'grRateKO', 'grRateWT', 'hasEffect', 'delRxns', 'fluxSolution'};
        else
            allowedSgd = {'grRatio', 'grRateKO', 'grRateWT', 'hasEffect', 'delRxns', 'fluxSolution'};
        end
        actualSgd = fieldnames(entry.singleGeneDeletion);
        extraSgd = setdiff(actualSgd, allowedSgd);
        if ~isempty(extraSgd)
            error("Unexpected field(s) in %s: %s. Allowed ones are: %s.", ...
                path, strjoin(extraSgd, ", "), strjoin(allowedSgd, ", "));
        end
        if isActive && isfield(entry.singleGeneDeletion, 'analysisId')
            if ~(ischar(entry.singleGeneDeletion.analysisId) || isstring(entry.singleGeneDeletion.analysisId))
                error("%s.singleGeneDeletion.analysisId must be a char or string.", path);
            end
        end
    end

    % Validate doubleGeneDeletion (struct with sub-fields + analysisId if active)
    if isfield(entry, 'doubleGeneDeletion')
        if ~isstruct(entry.doubleGeneDeletion)
            error("%s.doubleGeneDeletion must be a struct.", path);
        end
        if isActive
            allowedDgd = {'analysisId', 'grRatioDble', 'grRatioKO', 'grRateWT'};
        else
            allowedDgd = {'grRatioDble', 'grRatioKO', 'grRateWT'};
        end
        actualDgd = fieldnames(entry.doubleGeneDeletion);
        extraDgd = setdiff(actualDgd, allowedDgd);
        if ~isempty(extraDgd)
            error("Unexpected field(s) in %s: %s. Allowed ones are: %s.", ...
                path, strjoin(extraDgd, ", "), strjoin(allowedDgd, ", "));
        end
        if isActive && isfield(entry.doubleGeneDeletion, 'analysisId')
            if ~(ischar(entry.doubleGeneDeletion.analysisId) || isstring(entry.doubleGeneDeletion.analysisId))
                error("%s.doubleGeneDeletion.analysisId must be a char or string.", path);
            end
        end
    end

    % Validate sampling (nSamples may be empty — dimension checks will be skipped)
    if isfield(entry, 'sampling')
        if isempty(nSamples)
            warning("%s.sampling: could not determine number of samples from parameters table " + ...
                "(missing parameters table, no row with Parameter='options.nPointsReturned' " + ...
                "and Analysis='sampling', or unconvertible value). Dimension checks will be skipped.", ...
                path);
        end
        validateSamplingInAnalysis(entry.sampling, path + ".sampling", nbRxns, nSamples, isActive);
    end

    % Validate kld (struct with sub-fields + analysisId if active)
    if isfield(entry, 'kld')
        if ~isstruct(entry.kld)
            error("%s.kld must be a struct.", path);
        end
        if isActive
            allowedKld = {'analysisId', 'samplingSets', 'kldMatrix', 'pValueKld', 'fdr'};
        else
            allowedKld = {'samplingSets', 'kldMatrix', 'pValueKld', 'fdr'};
        end
        actualKld = fieldnames(entry.kld);
        extraKld = setdiff(actualKld, allowedKld);
        if ~isempty(extraKld)
            error("Unexpected field(s) in %s: %s. Allowed ones are: %s.", ...
                path, strjoin(extraKld, ", "), strjoin(allowedKld, ", "));
        end
        if isActive && isfield(entry.kld, 'analysisId')
            if ~(ischar(entry.kld.analysisId) || isstring(entry.kld.analysisId))
                error("%s.kld.analysisId must be a char or string.", path);
            end
        end
    end

end

function validateSamplingInAnalysis(s, path, nbRxns, nSamples, isActive)
% Validate the sampling struct within an analysis entry.
% If nSamples is empty, dimension checks are skipped with a warning.
% If isActive is true, 'analysisId' is allowed as a subfield.

    if ~isstruct(s)
        error("%s must be a struct.", path);
    end

    if isActive
        allowed = {'analysisId', 'modelSampling', 'samples', 'cycleFreeFlux'};
    else
        allowed = {'modelSampling', 'samples', 'cycleFreeFlux'};
    end
    actual = fieldnames(s);
    extra = setdiff(actual, allowed);
    if ~isempty(extra)
        error("Unexpected field(s) in %s: %s. Allowed ones are: %s.", ...
            path, strjoin(extra, ", "), strjoin(allowed, ", "));
    end

    % Validate analysisId (char, only in active)
    if isActive && isfield(s, 'analysisId')
        if ~(ischar(s.analysisId) || isstring(s.analysisId))
            error("%s.analysisId must be a char or string.", path);
        end
    end

    checkDims = ~isempty(nSamples);

    % Validate modelSampling (struct, no sub-field format specified)
    if isfield(s, 'modelSampling') && ~isstruct(s.modelSampling)
        error("%s.modelSampling must be a struct.", path);
    end

    % Validate samples (single, nbRxns x nSamples)
    if isfield(s, 'samples')
        if ~isa(s.samples, 'single')
            error("%s.samples must be single precision.", path);
        end
        if checkDims && ~isequal(size(s.samples), [nbRxns, nSamples])
            warning("%s.samples dimensions are %s but expected %dx%d (nbRxns x nbSamples). " + ...
                "Dimensions do not match the parameters table.", ...
                path, mat2str(size(s.samples)), nbRxns, nSamples);
        end
    end

    % Validate cycleFreeFlux
    if isfield(s, 'cycleFreeFlux')
        validateCycleFreeFluxInAnalysis(s.cycleFreeFlux, path + ".cycleFreeFlux", nbRxns, nSamples);
    end

end

function validateCycleFreeFluxInAnalysis(cff, path, nbRxns, nSamples)
% Validate the cycleFreeFlux struct within sampling.
% If nSamples is empty, dimension checks are skipped with a warning.

    if ~isstruct(cff)
        error("%s must be a struct.", path);
    end

    allowed = {'samplesLl', 'thermoFeas', 'sampleStatusAfterCorrection', ...
        'neededAttempts', 'looplessStatus'};
    actual = fieldnames(cff);
    extra = setdiff(actual, allowed);
    if ~isempty(extra)
        error("Unexpected field(s) in %s: %s. Allowed ones are: 'samplesLl', 'thermoFeas', " + ...
            "'sampleStatusAfterCorrection', 'neededAttempts', 'looplessStatus'.", ...
            path, strjoin(extra, ", "));
    end

    checkDims = ~isempty(nSamples);

    % Validate samplesLl (single, nbRxns x nSamples)
    if isfield(cff, 'samplesLl')
        if ~isa(cff.samplesLl, 'single')
            error("%s.samplesLl must be single precision.", path);
        end
        if checkDims && ~isequal(size(cff.samplesLl), [nbRxns, nSamples])
            warning("%s.samplesLl dimensions are %s but expected %dx%d (nbRxns x nbSamples). " + ...
                "Dimensions do not match the parameters table.", ...
                path, mat2str(size(cff.samplesLl)), nbRxns, nSamples);
        end
    end

    % Validate thermoFeas (uint8, nbRxns x nSamples)
    if isfield(cff, 'thermoFeas')
        if ~isa(cff.thermoFeas, 'uint8')
            error("%s.thermoFeas must be uint8.", path);
        end
        if checkDims && ~isequal(size(cff.thermoFeas), [nbRxns, nSamples])
            warning("%s.thermoFeas dimensions are %s but expected %dx%d (nbRxns x nbSamples). " + ...
                "Dimensions do not match the parameters table.", ...
                path, mat2str(size(cff.thermoFeas)), nbRxns, nSamples);
        end
    end

    % Validate sampleStatusAfterCorrection (uint8, nbRxns x nSamples)
    if isfield(cff, 'sampleStatusAfterCorrection')
        if ~isa(cff.sampleStatusAfterCorrection, 'uint8')
            error("%s.sampleStatusAfterCorrection must be uint8.", path);
        end
        if checkDims && ~isequal(size(cff.sampleStatusAfterCorrection), [nbRxns, nSamples])
            warning("%s.sampleStatusAfterCorrection dimensions are %s but expected %dx%d (nbRxns x nbSamples). " + ...
                "Dimensions do not match the parameters table.", ...
                path, mat2str(size(cff.sampleStatusAfterCorrection)), nbRxns, nSamples);
        end
    end

    % Validate neededAttempts (uint8, nbRxns x 1)
    if isfield(cff, 'neededAttempts')
        if ~isa(cff.neededAttempts, 'uint8')
            error("%s.neededAttempts must be uint8.", path);
        end
        if checkDims && ~isequal(size(cff.neededAttempts), [nbRxns, 1])
            warning("%s.neededAttempts dimensions are %s but expected %dx1 (nbRxns x 1). " + ...
                "Dimensions do not match the parameters table.", ...
                path, mat2str(size(cff.neededAttempts)), nbRxns);
        end
    end

    % Validate looplessStatus (uint8, nbRxns x nSamples)
    if isfield(cff, 'looplessStatus')
        if ~isa(cff.looplessStatus, 'uint8')
            error("%s.looplessStatus must be uint8.", path);
        end
        if checkDims && ~isequal(size(cff.looplessStatus), [nbRxns, nSamples])
            warning("%s.looplessStatus dimensions are %s but expected %dx%d (nbRxns x nbSamples). " + ...
                "Dimensions do not match the parameters table.", ...
                path, mat2str(size(cff.looplessStatus)), nbRxns, nSamples);
        end
    end

end

function nSamples = getNumberOfSamples(parameters)
% Extract the number of samples from the parameters table.
% Looks for the row where Parameter = 'options.nPointsReturned'
% and Analysis = 'sampling', and converts the Value (char) to a number.
% Returns [] if the row is not found or the value cannot be converted.

    nSamples = [];

    mask = (string(parameters.Parameter) == "options.nPointsReturned") & ...
           (string(parameters.Analysis) == "sampling");

    if any(mask)
        idx = find(mask, 1);
        val = parameters.Value{idx};
        nSamples = str2double(val);
        if isnan(nSamples)
            warning("Could not convert value '%s' to a number in parameters " + ...
                "(options.nPointsReturned/sampling). Dimension checks will be skipped.", val);
            nSamples = [];
        end
    end

end

%%
% project.models
% project.models.Name1
% project.models.Name1.model
% project.models.Name1.sampleMetadata
% project.models.Name1.discretizedData
% project.models.Name1.expressionData
% project.models.Name1.geneIds
% project.models.Name1.mappedDiscretizedRxns
% project.models.Name1.mappingDiscretizedRxnsAllSamples
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

% project.models.Name1.analysis (struct)
% project.models.Name1.analysis.analysis_id (struct)
% project.models.Name1.analysis.analysis_id.parameters (table with Parameter, Analysis and Value column)
% project.models.Name1.analysis.analysis_id.FBA (struct)
% project.models.Name1.analysis.analysis_id.FVA (struct)
% project.models.Name1.analysis.analysis_id.FVA.minMaxFluxes (table with minFlux, maxFlux, rows = nb of rxns)
% project.models.Name1.analysis.analysis_id.FVA.loopStatus (logical, dim nb of rxns*1)
% project.models.Name1.analysis.analysis_id.singleGeneDeletion (struct)
% project.models.Name1.analysis.analysis_id.singleGeneDeletion.grRatio
% project.models.Name1.analysis.analysis_id.singleGeneDeletion.grRateKO
% project.models.Name1.analysis.analysis_id.singleGeneDeletion.grRateWT
% project.models.Name1.analysis.analysis_id.singleGeneDeletion.hasEffect
% project.models.Name1.analysis.analysis_id.singleGeneDeletion.delRxns
% project.models.Name1.analysis.analysis_id.singleGeneDeletion.fluxSolution
% project.models.Name1.analysis.analysis_id.doubleGeneDeletion (struct)
% project.models.Name1.analysis.analysis_id.doubleGeneDeletion.grRatioDble
% project.models.Name1.analysis.analysis_id.doubleGeneDeletion.grRatioKO
% project.models.Name1.analysis.analysis_id.doubleGeneDeletion.grRateWT
% project.models.Name1.analysis.analysis_id.sampling (struct)
% project.models.Name1.analysis.analysis_id.sampling.modelSampling (struct)
% project.models.Name1.analysis.analysis_id.sampling.samples (single, dim nb of rxns * nb of samples)
% project.models.Name1.analysis.analysis_id.sampling.cycleFreeFlux (struct)
% project.models.Name1.analysis.analysis_id.sampling.cycleFreeFlux.samplesLl (single, dim nb of rxns * nb of samples)
% project.models.Name1.analysis.analysis_id.sampling.cycleFreeFlux.thermoFeas (uint8, dim nb of rxns * nb of samples)
% project.models.Name1.analysis.analysis_id.sampling.cycleFreeFlux.sampleStatusAfterCorrection (uint8, dim nb of rxns * nb of samples);
% project.models.Name1.analysis.analysis_id.sampling.cycleFreeFlux.neededAttempts (uint8, nb of rxns * 1);
% project.models.Name1.analysis.analysis_id.sampling.cycleFreeFlux.looplessStatus (uint8, dim nb of rxns * nb of samples);
% project.models.Name1.analysis.analysis_id.kld (struct)
% project.models.Name1.analysis.analysis_id.kld.samplingSets
% project.models.Name1.analysis.analysis_id.kld.kldMatrix
% project.models.Name1.analysis.analysis_id.kld.pValueKld
% project.models.Name1.analysis.analysis_id.kld.fdr

% project.models.Name1.analysis.active (struct)
% project.models.Name1.analysis.active.parameters (table with Parameter, Analysis and Value column)
% project.models.Name1.analysis.active.FBA (struct)
% project.models.Name1.analysis.active.FBA.analysisId
% project.models.Name1.analysis.active.FVA (struct) 
% project.models.Name1.analysis.active.FVA.minMaxFluxes (table with minFlux, maxFlux, rows = nb of rxns)
% project.models.Name1.analysis.active.FVA.loopStatus (logical, dim nb of rxns*1)
% project.models.Name1.analysis.active.FVA.analysisId (char)
% project.models.Name1.analysis.active.singleGeneDeletion (struct)
% project.models.Name1.analysis.active.singleGeneDeletion.analysisId
% project.models.Name1.analysis.active.singleGeneDeletion.grRatio
% project.models.Name1.analysis.active.singleGeneDeletion.grRateKO
% project.models.Name1.analysis.active.singleGeneDeletion.grRateWT
% project.models.Name1.analysis.active.singleGeneDeletion.hasEffect
% project.models.Name1.analysis.active.singleGeneDeletion.delRxns
% project.models.Name1.analysis.active.singleGeneDeletion.fluxSolution
% project.models.Name1.analysis.active.doubleGeneDeletion (struct)
% project.models.Name1.analysis.active.doubleGeneDeletion.analysisId
% project.models.Name1.analysis.active.doubleGeneDeletion.grRatioDble
% project.models.Name1.analysis.active.doubleGeneDeletion.grRatioKO
% project.models.Name1.analysis.active.doubleGeneDeletion.grRateWT
% project.models.Name1.analysis.active.sampling (struct)
% project.models.Name1.analysis.active.sampling.analysisId
% project.models.Name1.analysis.active.sampling.modelSampling (struct)
% project.models.Name1.analysis.active.sampling.samples (single, dim nb of rxns * nb of samples)
% project.models.Name1.analysis.active.sampling.cycleFreeFlux (struct)
% project.models.Name1.analysis.active.sampling.cycleFreeFlux.samplesLl (single, dim nb of rxns * nb of samples)
% project.models.Name1.analysis.active.sampling.cycleFreeFlux.thermoFeas (uint8, dim nb of rxns * nb of samples)
% project.models.Name1.analysis.active.sampling.cycleFreeFlux.sampleStatusAfterCorrection (uint8, dim nb of rxns * nb of samples);
% project.models.Name1.analysis.active.sampling.cycleFreeFlux.neededAttempts (uint8, nb of rxns * 1);
% project.models.Name1.analysis.active.sampling.cycleFreeFlux.looplessStatus (uint8, dim nb of rxns * nb of samples);
% project.models.Name1.analysis.active.kld (struct)
% project.models.Name1.analysis.active.kld.analysisId
% project.models.Name1.analysis.active.kld.samplingSets
% project.models.Name1.analysis.active.kld.kldMatrix
% project.models.Name1.analysis.active.kld.pValueKld
% project.models.Name1.analysis.active.kld.fdr

% project.comparisons (struct)
% project.comparisons.Name1_vs_Name2__date (struct)
% project.comparisons.Name1_vs_Name2_date.structuralAnalysisStatus (maybe not needed anymore if we use checkProjectFormat instead ?)
% project.comparisons.Name1_vs_Name2_date.modelNames (nb of models compared * 1, string)
% project.comparisons.Name1_vs_Name2_date.rxnMappingTable (table, nb of rxns in ref model * nb of models compared)
% project.comparisons.Name1_vs_Name2_date.referenceModel (name of ref model, string) (needs to be in project.models)
% project.comparisons.Name1_vs_Name2_date.orderedFBA. (double, nb of rxns in ref model * nb of models compared)
% project.comparisons.Name1_vs_Name2_date.orderedSamples (double, nb of rxn in ref model * cumulative nb of sampels over compared models)
% project.comparisons.Name1_vs_Name2_date.sampleModelLabels (string, dim 1*cumulative nb of sampels over compared models)
% project.comparisons.Name1_vs_Name2_date.plots (struct)

