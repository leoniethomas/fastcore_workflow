function [rxnID,producing,matchedAll] = getRxnIDs(project,referenceModel, pattern)
    % This function allow you to filter for reactions to visualize using a
    % pattern. You can give in the pattern for a subsystem/gene/rxns or
    % metabolite and the function will find the rxns that are associated
    % with the pattern given as an input. 
    % Additionally you can define multiple patterns to be matched. 
    % If you want the rxns to be associated to both patterns (for example a
    % rxns that produces/consumes lactate  in the Glycolysis) then you can
    % give in a pattern like this "^lac_.* & ^Glycolysis.*". If you want
    % to find rxns that are assocated with multiple patterns (for example
    % all reactions that are consuming/producing lactate and pyruvate) you
    % can write it like this: "^lac.* | ^pyr.*". 
    % USAGE:
    %
    %   [rxnID,producing,matched] = getRxnIDs(project,referenceModel, pattern)
    % INPUTS:
    %   project:                project object which is the output of
    %                           singleModelAnalysis and modelsComparison
    %   referenceModel:         the reference model according to which the
    %                           idx should be obtained
    %   pattern:                the regex pattern that is matched to find
    %                           genes/rxns/mets participating in a rxns
    % OUTPUTS: 
    %   rxnID:                  idx of the models in the reference model
    %   producing:              logical vector indicating if the mets
    %                           specified in the pattern are produced or consumed in the reaction
    %                           (according to stochiometry)
    %   matched:                returning the actuall strings the pattern
    %                           were matched to (so the name of the
    %                           rxns/genes/mets/subSystems detected by
    %                           givin the pattern)
    arguments
        project (1,1) struct
        referenceModel (1,1) string
        pattern (1,:) string
    end

    model = project.models.(referenceModel).model;
    dict = project.models.(referenceModel).settings.dico;
    rxnID = {};
    producing = {};
    matchedAll = {};

    for n= 1:numel(pattern)
        
        singlePat = pattern(n);

        assert(~(contains(singlePat, "|") & contains(singlePat, "&")),...
               'The pattern can not contain both & + |, this function is only written for either or !')
            
        if contains(singlePat, "|")
            patCond =  strtrim(strsplit(singlePat, "|"));
        elseif contains(singlePat, "&")
            patCond =  strtrim(strsplit(singlePat, "&"));  
        else
            patCond = strtrim(singlePat);
        end
        
        % filter out all fields that can not be matched to pattern 
        fieldsToCheck = structfun(@(x) (ischar(x) || isstring(x) || iscell(x) ) ...
                                        && ((size(x,2) == 1 && size(x,1) == size(model.mets,1) || size(x,1) == size(model.rxns,1)|| size(x,1) == size(model.genes,1))),...
                                  model);
        namesToCheck = fieldnames(model);
        namesToNotCheck = namesToCheck(~fieldsToCheck);
        slotsCheck = rmfield(model,namesToNotCheck);
    
        [rxnsIDsAll,producingMetAll,matched] = cellfun(@(pat) findPatternsInStruct(pat,slotsCheck,model,dict),...
                                               patCond, 'UniformOutput', false);
    
        if contains(singlePat, "&")
            commonRxns = rxnsIDsAll{1};
            commonProd = producingMetAll{1};
            for i = 2:numel(rxnsIDsAll)
                [commonRxns, ia] = intersect(commonRxns, rxnsIDsAll{i});
                commonProd = commonProd(ia);
            end
            resultRxns = commonRxns;
            resultProd = commonProd;
            
        elseif contains(singlePat, "|")
            allRxns = vertcat(rxnsIDsAll{:});
            allProd = vertcat(producingMetAll{:});
            [resultRxns, ia] = unique(allRxns, 'stable');
            resultProd = allProd;
        else
            resultRxns = rxnsIDsAll{:,:};
            resultProd = producingMetAll{:,:};
            
        end
        
        rxnID{n} = resultRxns;
        producing{n} = resultProd;
        matchedAll{n} = string(vertcat(matched{:}));
        assert(~isempty(matchedAll{n}), 'No reactions were found! Check your pattern for spelling mistakes!')
    end
end

function [rxnIDs, producingMet,matched] = findPatternsInStruct(pattern,slotsCheck,model,dico)
    SMatrix = full(model.S);
    patMatch = structfun(@(slot) regexp(string(slot), pattern, 'once'),...
                         slotsCheck, 'UniformOutput',false);
    posMatch = structfun(@(x) find(~cellfun(@isempty, x)), ...
                     patMatch, 'UniformOutput', false);

    if all(structfun(@length, posMatch) == 0)
        % this means that the string given is a gene ? 
        % can we find a match in the dico ? 
        matches = findInDico(dico, pattern);
        colname = getDicoColumnWithGeneIDsInModel(model, dico);
        if length(matches.rowIdx) ~= 0
            genesPattern = matches.rows.(colname);
            model.genes = regexprep(model.genes, '\.\d+$', '');
            rxnNames = findRxnsFromGenes(model,genesPattern);
            matched = string(fieldnames(rxnNames));
            rxnNames = struct2cell(structfun(@(x) x(:,1), rxnNames, 'UniformOutput', false));
            allNames = string(unique(vertcat(rxnNames{:})));
            [~,rxnIDs] = ismember(allNames,model.rxns);
            producingMet = repmat(0, 1,length(rxnIDs));
        else
            matched = [];
            rxnIDs = [];
            producingMet = [];
        end
    
    else
        [~,mostMatchSlot] = max(structfun(@(x) size(x,1),posMatch));
        metdimStatus = structfun(@(slot) size(slot,1) == size(slotsCheck.mets,1),slotsCheck);
        patternDefMets = metdimStatus(mostMatchSlot);
        fieldName = fieldnames(posMatch);
        fieldName = string(fieldName(mostMatchSlot));
    
        idxMatched = posMatch.(fieldName);
        matched = slotsCheck.(fieldName);
        matched = matched(idxMatched);
    
        % now we have the IDx of the matches, if the matched slot is in the
        % dimension of metabolites, the rxns have to be found from the
        % metabolite 
    
        if patternDefMets
            % find metablite names from the idx of the matches 
            metNames = slotsCheck.mets(idxMatched);
            rxnNames = findRxnsFromMets(model, metNames);
            [~,rxnIDs] = ismember(rxnNames, slotsCheck.rxns);
            producingMet = any(SMatrix(idxMatched,rxnIDs) == 1);
        else
            rxnIDs = idxMatched;
            producingMet = repmat(0, 1,length(rxnIDs));
        end

    end
    producingMet = producingMet';
end


function matches = findInDico(dico, query)
% FINDINDICO  Find a query string in any column of the dico table.
%   matches = findInDico(dico, query)
%
%   Returns a struct with:
%     matches.rows      - subtable of matching rows
%     matches.columns   - cell array of column names where match was found
%     matches.rowIdx    - row indices in original table

    colNames = dico.Properties.VariableNames; % {'ENTREZ','HGNC','ENSG','SYMBOL'}
    
    nCols = numel(colNames);
    matchMask = false(height(dico), nCols);
    
    for c = 1:nCols
        col = dico.(colNames{c});
        % Each column is a cell array of char vectors
        matchMask(:, c) = ~cellfun(@isempty, regexp(col, query, 'once'));
    end
    
    % Rows that matched in at least one column
    rowHits  = any(matchMask, 2);
    colHits  = colNames(any(matchMask, 1));
    
    matches.rows    = dico(rowHits, :);
    matches.columns = colHits;
    matches.rowIdx  = find(rowHits);
end

function colname = getDicoColumnWithGeneIDsInModel(model, dico)
        queryRaw = string(model.genes); % your vector

        queryStripped = regexprep(queryRaw, '\.\d+$', '');
        
        colNames = dico.Properties.VariableNames;
        overlapCount = zeros(1, numel(colNames));
        
        for c = 1:numel(colNames)
            col = dico.(colNames{c});
            overlapCount(c) = max(numel(intersect(queryRaw, col)), ...
                                 numel(intersect(queryStripped, col)));
        end
        
        [~, bestIdx] = max(overlapCount);
        colname = colNames{bestIdx};

end