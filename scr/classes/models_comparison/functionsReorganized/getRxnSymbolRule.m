function rxnGeneRule = getRxnSymbolRule(model, rxnName)
    % --- Get gene association for reaction ---
    % Replace gene specifiers by gene names if geneNames column exists,
    % otherwise keep original IDs from geneIdsInModel

    % --- Locate reaction ---
    rxnIdx = find(strcmp(model.model.rxns, rxnName), 1);
    assert(~isempty(rxnIdx), 'Reaction not found in the model.');

    % --- Retrieve gene association string ---
    geneAssociation = model.model.grRules{rxnIdx};
    geneTable = model.settings.dico;

    % --- Extract gene IDs with isoform suffix (e.g. 1234.1) ---
    idsWithPostfix = regexp(geneAssociation, '\d+\.\d+', 'match');

    % No gene IDs in this reaction: return as-is
    if isempty(idsWithPostfix)
        rxnGeneRule = geneAssociation;
        return;
    end

    % --- Check if geneNames column exists ---
    useGeneNames = ismember('geneNames', geneTable.Properties.VariableNames);

    % No geneNames column: keep original IDs unchanged
    if ~useGeneNames
        rxnGeneRule = geneAssociation;
        return;
    end

    % --- Map to geneIdsInModel ---
    [~, idx] = ismember(string(idsWithPostfix), string(geneTable.geneIdsInModel));

    % --- Build symbols: geneNames if available, else geneIdsInModel ---
    geneNamesCol = string(geneTable.geneNames);
    geneIdsCol = string(geneTable.geneIdsInModel);

    valid = idx > 0;
    symbols = repmat(string(""), numel(idsWithPostfix), 1);

    % Fallback: use geneIdsInModel for all valid entries
    symbols(valid) = geneIdsCol(idx(valid));

    % Override with geneNames where non-empty and not missing
    mappedNames = geneNamesCol(idx(valid));
    namesValid = ~ismissing(mappedNames) & strlength(mappedNames) > 0;
    validIdx = find(valid);
    symbols(validIdx(namesValid)) = mappedNames(namesValid);

    % --- Substitute IDs by symbols in the association string ---
    geneAssocOut = geneAssociation;

    for k = 1:numel(idsWithPostfix)
        if strlength(symbols(k)) > 0
            geneAssocOut = regexprep(geneAssocOut, ...
                ['\<' idsWithPostfix{k} '\>'], char(symbols(k)));
        else
            warning('ID %s not found in geneIdsInModel; keeping original ID.', ...
                idsWithPostfix{k});
        end
    end

    rxnGeneRule = geneAssocOut;

end