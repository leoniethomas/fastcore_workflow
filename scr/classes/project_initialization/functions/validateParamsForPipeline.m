function params = validateParamsForPipeline(params)
% Validates a paramsForModel structure.
% Takes a structure and validates that all required and optional fields
% conform to the expected types, formats, and constraints for the pipeline.

    arguments
        params struct
    end

    % Required fields
    
    if ~isfield(params, 'modelName')
        error("params.modelName is required.");
    end
    mustBeText(params.modelName);

    if ~isfield(params, 'contextSpecificModel')
        error("params.contextSpecificModel is required.");
    end
    if ~(isscalar(params.contextSpecificModel) && isa(params.contextSpecificModel, 'struct'))
        error("params.contextSpecificModel must be a scalar struct.");
    end
    mustBeCobraModel(params.contextSpecificModel);
    
    % If provided, validation of optional fields

    if isFieldNonEmpty(params, 'dico') && ~isa(params.dico, 'table')
        error("params.dico must be a table.");
    end

    if isFieldNonEmpty(params, 'sampleMetadata') && ~isa(params.sampleMetadata, 'table')
        error("params.sampleMetadata must be a table.");
    end

    if isFieldNonEmpty(params, 'sampleLabeling')
        mustBeText(params.sampleLabeling);
    end

    if isFieldNonEmpty(params, 'mediumComposition') && ~isa(params.mediumComposition, 'table')
        error("params.mediumComposition must be a table.");
    end

    if isFieldNonEmpty(params, 'objFunction')
        mustBeText(params.objFunction);
    end

    if isFieldNonEmpty(params, 'discretizedData')
        if ~isa(params.discretizedData, 'double') || ndims(params.discretizedData) ~= 2
            error("params.discretizedData must be a 2-D double matrix.");
        end
    end

    if isFieldNonEmpty(params, 'expressionData')
        if ~isa(params.expressionData, 'double') || ndims(params.expressionData) ~= 2
            error("params.expressionData must be a 2-D double matrix.");
        end
    end

    if isFieldNonEmpty(params, 'geneIds')
        % Accepte char / cellstr et normalise en string
        if isa(params.geneIds, 'char') || ...
           (isa(params.geneIds, 'cell') && all(cellfun(@ischar, params.geneIds(:))))
            params.geneIds = string(params.geneIds);
        end
        if ~isa(params.geneIds, 'string')
            error("params.geneIds must be a string array (char/cellstr accepted and converted).");
        end
    end

    if isFieldNonEmpty(params, 'consensusProportion')
        mustBeNonnegative(params.consensusProportion);
    end

    if isFieldNonEmpty(params, 'referenceModel')
        mustBeText(params.referenceModel);
    end

    if isFieldNonEmpty(params, 'optionalSettings')
        if ~(isscalar(params.optionalSettings) && isa(params.optionalSettings, 'struct'))
            error("params.optionalSettings must be a scalar struct.");
        end
        validateOptionalSettings(params.optionalSettings);
    end

    if isFieldNonEmpty(params, 'coreReactions')
        if ~isa(params.coreReactions, 'cell')
            error("params.coreReactions must be a cell array.");
        end
        mustBeVectorOrEmpty(params.coreReactions);
    end

    if isFieldNonEmpty(params, 'manuallySetBoundaries')
        if ~(isscalar(params.manuallySetBoundaries) && isa(params.manuallySetBoundaries, 'struct'))
            error("params.manuallySetBoundaries must be a scalar struct.");
        end
        validateManuallySetBoundaries(params.manuallySetBoundaries);
    end

    if isFieldNonEmpty(params, 'mapping')
        if ~isa(params.mapping, 'double')
            error("params.mapping must be a double matrix (sparse or dense).");
        end
    end

    % Cross-fields validation

    hasSampleMetadata = isFieldNonEmpty(params, 'sampleMetadata');
    hasSampleLabeling = isFieldNonEmpty(params, 'sampleLabeling') && ...
        any(strlength(string(params.sampleLabeling)) > 0);

    if hasSampleMetadata ~= hasSampleLabeling
        error("sampleMetadata and sampleLabeling must be provided together or both omitted.");
    end

    if hasSampleMetadata && ~ismember(string(params.sampleLabeling), ...
                                      string(params.sampleMetadata.Properties.VariableNames))
        error("'%s' must be a column in sampleMetadata.", params.sampleLabeling);
    end

    hasDiscretized = isFieldNonEmpty(params, 'discretizedData');
    hasExpression  = isFieldNonEmpty(params, 'expressionData');
    hasData        = hasDiscretized || hasExpression;
    hasGeneIds     = isFieldNonEmpty(params, 'geneIds') && ...
                     any(strlength(string(params.geneIds)) > 0);
    hasDico        = isFieldNonEmpty(params, 'dico');

    if hasData && ~hasGeneIds
        error("geneIds is required when expressionData or discretizedData is provided.");
    end

    if hasData && ~hasDico
        error("dico is required when expressionData or discretizedData is provided.");
    end

    if hasGeneIds && hasData
        if hasDiscretized && size(params.geneIds, 1) ~= size(params.discretizedData, 1)
            error("geneIds must have the same number of rows as discretizedData (%d rows), got %d rows.", ...
                size(params.discretizedData, 1), size(params.geneIds, 1));
        end
        if hasExpression && size(params.geneIds, 1) ~= size(params.expressionData, 1)
            error("geneIds must have the same number of rows as expressionData (%d rows), got %d rows.", ...
                size(params.expressionData, 1), size(params.geneIds, 1));
        end
    end

    if hasDiscretized && hasExpression && ...
       ~isequal(size(params.discretizedData), size(params.expressionData))
        error("discretizedData (%s) and expressionData (%s) must have the same dimensions (rows = genes, columns = samples).", ...
            mat2str(size(params.discretizedData)), mat2str(size(params.expressionData)));
    end

    if hasDico
        requiredDicoCols = {'geneIdsInModel', 'geneIdsInData'};
        missingCols = setdiff(requiredDicoCols, string(params.dico.Properties.VariableNames));
        if ~isempty(missingCols)
            error("dico must contain columns: %s. Missing: %s.", ...
                strjoin(requiredDicoCols, ", "), strjoin(missingCols, ", "));
        end

        modelGenes     = regexprep(string(params.contextSpecificModel.genes), "\.[0-9]+$", "");
        dicoModelGenes = regexprep(string(params.dico.geneIdsInModel), "\.[0-9]+$", "");
        nModelGenes    = numel(modelGenes);
        nOverlapModel  = sum(ismember(dicoModelGenes, modelGenes));
        pctModel       = (nOverlapModel / nModelGenes) * 100;

        if pctModel < 5
            error("dico.geneIdsInModel covers only %.1f%% (%d/%d) of contextSpecificModel genes. Minimum 5%% required.", ...
                pctModel, nOverlapModel, nModelGenes);
        end
        warning("dico.geneIdsInModel covers %.1f%% (%d/%d) of contextSpecificModel genes.", ...
            pctModel, nOverlapModel, nModelGenes);

        if hasData
            dataGeneIds   = regexprep(string(params.geneIds), "\.[0-9]+$", "");
            dicoDataGenes = regexprep(string(params.dico.geneIdsInData), "\.[0-9]+$", "");
            nDataGenes    = numel(dataGeneIds);
            nOverlapData  = sum(ismember(dicoDataGenes, dataGeneIds));
            pctData       = (nOverlapData / nDataGenes) * 100;

            if pctData < 5
                error("dico.geneIdsInData covers only %.1f%% (%d/%d) of geneIds. Minimum 5%% required.", ...
                    pctData, nOverlapData, nDataGenes);
            end
            warning("dico.geneIdsInData covers %.1f%% (%d/%d) of geneIds.", ...
                pctData, nOverlapData, nDataGenes);
        end
    end
end

function tf = isFieldNonEmpty(params, name)
% Return true if field exists and is not empty
    tf = isfield(params, name) && ~isempty(params.(name));
end