function modelShortcut = formatParamsForModel(paramsForModel)
% Reformats pipeline parameters into a model struct.
% Takes a paramsForModel structure (as validated by validateParamsForPipeline)
% and reorganizes its fields into the nested structure format expected by the
% pipeline, separating model data, medium settings, and script parameters.

modelShortcut = struct();
modelShortcut.model = paramsForModel.contextSpecificModel;

% Optional fields
if length(fieldnames(paramsForModel)) > 2
    modelShortcut.settings = struct();
    modelShortcut.settings.medium = struct();
    modelShortcut.settings.scriptParameters = struct();
    modelShortcut.settings.medium.manuallySetBoundaries = struct();

    % --- dico ---
    if isfield(paramsForModel, 'dico')
        dico = paramsForModel.dico;
        % Reorder dico so geneIdsInModel follows the order of model.genes
        % Genes in model but not in dico -> "NA" rows
        % Genes in dico but not in model -> dropped

        % Strip version suffixes (e.g. ".1") for matching
        modelGenes = regexprep(string(modelShortcut.model.genes), "\.[0-9]+$", "");
        dicoGenes = regexprep(string(dico.geneIdsInModel), "\.[0-9]+$", "");

        [isFound, idxDico] = ismember(modelGenes, dicoGenes);

        % Build aligned table: one row per model gene, NA where not found in dico
        dicoAligned = array2table(repmat("NA", numel(modelGenes), width(dico)), ...
            'VariableNames', dico.Properties.VariableNames);

        dicoAligned(isFound, :) = varfun(@string, dico(idxDico(isFound), :));

        modelShortcut.settings.dico = dicoAligned;
    else
        modelShortcut.settings.dico = table(modelShortcut.model.genes, 'VariableNames', 'geneIdsInModel');
    end

    % --- sampleMetadata ---
    if isfield(paramsForModel, 'sampleMetadata')
        modelShortcut.sampleMetadata = paramsForModel.sampleMetadata;
    end

    % --- discretizedData ---
    if isfield(paramsForModel, 'discretizedData')
        discretizedDataTable = table(string(paramsForModel.geneIds), ...
            paramsForModel.discretizedData, ...
            'VariableNames', {'geneIds', 'value'});

        discretizedDataAligned = reorderDataToDico(discretizedDataTable, ...
            modelShortcut.settings.dico, 'filter');

        modelShortcut.discretizedData = discretizedDataAligned;
    end

    % --- expressionData ---
    if isfield(paramsForModel, 'expressionData')
        expressionDataTable = table(string(paramsForModel.geneIds), ...
            paramsForModel.expressionData, ...
            'VariableNames', {'geneIds', 'expression'});

        expressionDataAligned = reorderDataToDico(expressionDataTable, ...
            modelShortcut.settings.dico, 'keep');

        modelShortcut.expressionData = expressionDataAligned;
    end

    % --- coreReactions ---
    if isfield(paramsForModel, 'coreReactions')
        modelShortcut.coreReactions = paramsForModel.coreReactions;
    end

    % --- mediumComposition ---
    if isfield(paramsForModel, 'mediumComposition')
        modelShortcut.settings.medium.mediumComposition = paramsForModel.mediumComposition;
    end

    % --- manuallySetBoundaries ---
    if isfield(paramsForModel, 'manuallySetBoundaries')
        modelShortcut.settings.medium.manuallySetBoundaries = paramsForModel.manuallySetBoundaries;
        if isfield(paramsForModel, 'optionalSettings')
            if isfield(paramsForModel.optionalSettings, 'notMediumConstrained')
                modelShortcut.settings.medium.manuallySetBoundaries.unconstrainedImports = unique([paramsForModel.manuallySetBoundaries.unconstrainedImports, ...
                    paramsForModel.optionalSettings.notMediumConstrained]);
            end
        end
    end

    % --- consensusProportion ---
    if isfield(paramsForModel, 'consensusProportion')
        modelShortcut.settings.scriptParameters.consensusProportion = paramsForModel.consensusProportion;
    end

    % --- optionalSettings ---
    if isfield(paramsForModel, 'optionalSettings')
        modelShortcut.settings.optionalSettings = paramsForModel.optionalSettings;
        if isfield(paramsForModel.optionalSettings, 'notMediumConstrained') && ~isfield(paramsForModel.manuallySetBoundaries, 'unconstrainedImports')
            modelShortcut.settings.medium.manuallySetBoundaries.unconstrainedImports = paramsForModel.optionalSettings.notMediumConstrained;
        end
    end

    % --- objFunction ---
    if isfield(paramsForModel, 'objFunction')
        modelShortcut.settings.objFunction = paramsForModel.objFunction;
    end

    % --- referenceModel ---
    if isfield(paramsForModel, 'referenceModel')
        modelShortcut.settings.referenceModel = paramsForModel.referenceModel;
    end

    % --- mapping ---
    if isfield(paramsForModel, 'mapping')
        mapping = paramsForModel.mapping;
        modelShortcut.settings.mapping = mapping; % previously mappedDiscRxnsSample
        mapping = full(mapping);
        if isfield(modelShortcut.settings.scriptParameters, "consensusProportion")
            numberOfSamples = size(mapping, 2);
            % definition of initialCore reactions
            modelShortcut.mappedDiscretizedRxns = sum(mapping == 1, 2) >= (modelShortcut.settings.scriptParameters.consensusProportion * numberOfSamples);
            % definition of the notExpressed genes
            notExpressed = find(sum(mapping == -1, 2) >= (modelShortcut.settings.scriptParameters.consensusProportion * numberOfSamples));
            modelShortcut.mappedDiscretizedRxns = int32(modelShortcut.mappedDiscretizedRxns);
            modelShortcut.mappedDiscretizedRxns(notExpressed) = -1;
            modelShortcut.mappedDiscretizedRxns = int8(modelShortcut.mappedDiscretizedRxns); % needs less storage
        else
            modelShortcut.mappedDiscretizedRxns = int8(repmat(-13, size(mapping, 1), 1)); % previously mappedDiscRxns
        end
    end

    % --- sampleLabeling ---
    if isfield(paramsForModel, 'sampleLabeling')
        modelShortcut.settings.scriptParameters.sampleLabeling = paramsForModel.sampleLabeling;
    end

    % Remove empty fields (except in .model fields)
    modelShortcut = removeRecursivelyEmptyFields(modelShortcut, {'model'});
end

end


function dataAligned = reorderDataToDico(dataTable, dico, mode)
% Reorder a data table to match dico gene order.
%   mode = 'filter': output has height(dico) rows, NaN for missing genes (discretizedData)
%   mode = 'keep':   output keeps all rows, dico-matched first then extras (expressionData)
%   Adds geneNames column from dico if available.

    dicoGenesInData = regexprep(string(dico.geneIdsInData), "\.[0-9]+$", "");
    dataGenes = regexprep(string(dataTable.geneIds), "\.[0-9]+$", "");

    if strcmp(mode, 'filter')
        % Filter mode: output size = dico, NaN for missing 
        [isPresent, idx] = ismember(dicoGenesInData, dataGenes);

        nCols = size(dataTable.value, 2);
        outValues = NaN(numel(dicoGenesInData), nCols);
        outValues(isPresent, :) = dataTable.value(idx(isPresent), :);

        dataAligned = table(dicoGenesInData, outValues, ...
            'VariableNames', ["geneIds", "value"]);

    else % 'keep'
        % Keep mode: reorder only, no rows added or removed
        [inDico, dicoIdx] = ismember(dataGenes, dicoGenesInData);

        expInDicoRows = find(inDico);
        [~, sortOrder] = sort(dicoIdx(expInDicoRows));
        orderedInDico = expInDicoRows(sortOrder);

        expNotInDicoRows = find(~inDico);
        finalOrder = [orderedInDico; expNotInDicoRows];

        dataAligned = dataTable(finalOrder, :);
    end

    % Add geneNames column from dico if available
    if ismember('geneNames', dico.Properties.VariableNames)
        alignedGenes = regexprep(string(dataAligned.geneIds), "\.[0-9]+$", "");
        [~, idxDico] = ismember(alignedGenes, dicoGenesInData);

        geneNamesCol = repmat("NA", height(dataAligned), 1);
        found = idxDico > 0;
        geneNamesCol(found) = string(dico.geneNames(idxDico(found)));

        dataAligned = addvars(dataAligned, geneNamesCol, ...
            'Before', 1, 'NewVariableNames', 'geneNames');
    end
end