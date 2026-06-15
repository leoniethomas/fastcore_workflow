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

end