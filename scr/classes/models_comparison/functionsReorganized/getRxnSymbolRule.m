function rxn_gene_rule = getRxnSymbolRule(model,rxnName)
        % --- Get gene association for reaction ---
        % repace the gene specifiers by the gene symbol
        
        rxnIdx = find(strcmp(model.model.rxns, rxnName), 1);
        assert(~isempty(rxnIdx), 'Reaction not found in the model.');
        
        geneAssociation = model.model.grRules{rxnIdx};
        geneTable       = model.settings.dico;
        
        % --- Extract numeric gene IDs ---
        ids_with_postfix = regexp(geneAssociation, '\d+\.\d+', 'match');
        ids_numeric      = string(regexp(ids_with_postfix, '^\d+', 'match'));
        
        % --- Find best-matching ID column ---
        cols = {string(geneTable.ENTREZ), string(geneTable.HGNC), ...
                string(geneTable.ENSG),   string(geneTable.SYMBOL)};
        
        [~, idxMax] = max(cellfun(@(c) sum(ismember(ids_numeric, c)), cols));
        
        % --- Map to symbols ---
        [~, idx] = ismember(ids_numeric, geneTable{:, idxMax});
    
        % Check that idx is valid before using it
        if isempty(idx) || any(idx <= 0)
            symbols = string([]);   % or assign some default value, e.g., ""
        else
            symbols = string(geneTable.SYMBOL(idx));
        end

        %symbols  = string(geneTable.SYMBOL(idx));
        
        geneAssocOut = geneAssociation;  % start with original string

        for k = 1:numel(ids_with_postfix)
            if k <= numel(symbols) && ~isempty(symbols(k))
                geneAssocOut = regexprep(geneAssocOut, ['\<' ids_with_postfix{k} '\>'], symbols(k));
            else
                % optional: assign empty string or leave original ID
                % geneAssocOut = regexprep(geneAssocOut, ['\<' ids_with_postfix{k} '\>'], '');
                warning('No symbol for ID %s; skipping replacement.', ids_with_postfix{k});
            end
        end
        
        rxn_gene_rule = geneAssocOut;

end
