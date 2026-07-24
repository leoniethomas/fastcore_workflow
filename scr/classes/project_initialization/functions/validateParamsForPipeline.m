function params = validateParamsForPipeline(params)
% Validates a paramsForModel structure.
% Takes a structure and validates that all required and optional fields
% conform to the expected types, formats, and constraints for the pipeline.

    arguments
        params.modelName {mustBeText}
        params.contextSpecificModel (1,1) struct {mustBeCobraModel}
        
        params.dico table = []
        params.sampleMetadata table = table()
        params.sampleLabeling {mustBeText} = ''
        params.mediumComposition table = table()
        params.objFunction {mustBeText} = ''
        params.discretizedData (:,:) double = []
        params.expressionData table = table()
        params.geneIds string = ''
        params.consensusProportion = []
        params.referenceModel {mustBeText} = ''
        params.optionalSettings (1,1) struct {validateOptionalSettings} = struct()
        params.mapping = []
        params.coreReactions cell {mustBeVectorOrEmpty} = {}
        params.manuallySetBoundaries (1,1) struct {validateManuallySetBoundaries} = struct()
    end
    
    % Cross-field validation
    % Either both empty or both provided
    if isempty(params.sampleMetadata) ~= (strlength(params.sampleLabeling) == 0)
        error("sampleMetadata and sampleLabeling must be provided together or both omitted.");
    end
    
    % When provided, sampleLabeling must be a column of sampleMetadata
    if ~isempty(params.sampleMetadata)
        if ~ismember(string(params.sampleLabeling), string(params.sampleMetadata.Properties.VariableNames))
            error("'%s' must be a column in sampleMetadata.", params.sampleLabeling);
        end
    end

    if ~isempty(params.consensusProportion)
        mustBeNonnegative(params.consensusProportion);
    end
    
    % geneIds: required when expressionData or discretizedData is provided
    hasData = ~isempty(params.discretizedData) || ~isempty(params.expressionData);
    hasGeneIds = isfield(params, 'geneIds') && ~isempty(params.geneIds);
    hasDico = isfield(params, 'dico') && ~isempty(params.dico);
    
    if hasData && ~hasGeneIds
        error("geneIds is required when expressionData or discretizedData is provided.");
    end
    
    if hasData && ~hasDico
        error("dico is required when expressionData or discretizedData is provided.");
    end
    
    if hasGeneIds && hasData
        if ~isempty(params.discretizedData) && size(params.geneIds, 1) ~= size(params.discretizedData, 1)
            error("geneIds must have the same number of rows as discretizedData (%d rows), got %d rows.", ...
                size(params.discretizedData, 1), size(params.geneIds, 1));
        end
        if ~isempty(params.expressionData) && size(params.geneIds, 1) ~= size(params.expressionData, 1)
            error("geneIds must have the same number of rows as expressionData (%d rows), got %d rows.", ...
                size(params.expressionData, 1), size(params.geneIds, 1));
        end
    end
    
    % expressionData and discretizedData must have identical dimensions when both provided
    if ~isempty(params.discretizedData) && ~isempty(params.expressionData)
        if ~isequal(size(params.discretizedData), size(params.expressionData))
            error("discretizedData (%s) and expressionData (%s) must have the same dimensions (rows = genes, columns = samples).", ...
                mat2str(size(params.discretizedData)), mat2str(size(params.expressionData)));
        end
    end
    
    if isfield(params, 'dico')
        requiredDicoCols = {'geneIdsInModel', 'geneIdsInData'};
        missingCols = setdiff(requiredDicoCols, string(params.dico.Properties.VariableNames));
        
        if ~isempty(missingCols)
            error("dico must contain columns: %s. Missing: %s.", ...
                strjoin(requiredDicoCols, ", "), strjoin(missingCols, ", "));
        end
        
        % geneIdsInModel must cover >= 5% of contextSpecificModel genes
        modelGenes = string(params.contextSpecificModel.genes);
        dicoModelGenes = string(params.dico.geneIdsInModel);
        nModelGenes = numel(modelGenes);
        nOverlapModel = sum(ismember(dicoModelGenes, modelGenes));
        pctModel = (nOverlapModel / nModelGenes) * 100;
        
        if pctModel < 5
            error("dico.geneIdsInModel covers only %.1f%% (%d/%d) of contextSpecificModel genes. Minimum 5%% required.", ...
                pctModel, nOverlapModel, nModelGenes);
        end
        warning("dico.geneIdsInModel covers %.1f%% (%d/%d) of contextSpecificModel genes.", ...
            pctModel, nOverlapModel, nModelGenes);
        
        % geneIdsInData must cover >= 5% of geneIds when data is provided
        if hasData
            dataGeneIds = string(params.geneIds);
            dicoDataGenes = string(params.dico.geneIdsInData);
            nDataGenes = numel(dataGeneIds);
            nOverlapData = sum(ismember(dicoDataGenes, dataGeneIds));
            pctData = (nOverlapData / nDataGenes) * 100;
            
            if pctData < 5
                error("dico.geneIdsInData covers only %.1f%% (%d/%d) of geneIds. Minimum 5%% required.", ...
                    pctData, nOverlapData, nDataGenes);
            end
            warning("dico.geneIdsInData covers %.1f%% (%d/%d) of geneIds.", ...
                pctData, nOverlapData, nDataGenes);
        end
    end

end