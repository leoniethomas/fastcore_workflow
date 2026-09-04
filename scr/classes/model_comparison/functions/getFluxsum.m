function fluxsumCell = getFluxsum(project,comparison_name,metIdx,rxnIdx,slot,fluxSummedUp)
    arguments
        project
        comparison_name
        metIdx =[]
        rxnIdx =[]
        slot  (1,1) string {mustBeMember(slot,["orderedFba", "orderedSamples"])} ="orderedSamples" 
        fluxSummedUp {mustBeMember(fluxSummedUp, ["incoming","outgoing","reactions"])} ="incoming" 
    end
    
    reference = project.comparisons.(comparison_name).referenceModel;
    if isempty(rxnIdx)
        rxnIdx = {find(ones(length(project.models.(reference).model.rxns),1))};
    end
    
    if isempty(metIdx)
        metIdx = find(ones(length(project.models.(reference).model.mets),1));
    end
    
    
    samples = project.comparisons.(comparison_name).(slot);
    fullS = project.models.(reference).model.S;
    
    fluxsumCell = cellfun(@(x) getModelFluxsum(fullS(:,x), samples(x,:),fluxSummedUp), ...
                           rxnIdx,'UniformOutput', false);
    if fluxSummedUp ~= "reactions"
        % this is only needed in case the fluxsum is computed over all the
        % metabolites in the given reactions, as an additional filtering
        % step, in case of "reactions" only one value is returend per set
        % for the getModelFluxsum
        fluxsumCell = cellfun(@(f) f(metIdx,:), ...
                               fluxsumCell, 'UniformOutput', false);
    end

    

end
