function fluxsum_cell = get_fluxsum(project,comparison_name,met_idx,rxn_idx,slot,flux_summed_up)
    arguments
        project
        comparison_name
        met_idx =[]
        rxn_idx =[]
        slot  (1,1) string {mustBeMember(slot,["ordered_fba", "ordered_samples"])} ="ordered_samples" 
        flux_summed_up {mustBeMember(flux_summed_up, ["incoming","outgoing","reactions"])} ="incoming" 
    end
    
    reference = project.comparisons.(comparison_name).reference_model;
    if isempty(rxn_idx)
        rxn_idx = {find(ones(length(project.models.(reference).model.rxns),1))};
    end
    
    if isempty(met_idx)
        met_idx = find(ones(length(project.models.(reference).model.mets),1));
    end
    
    
    samples = project.comparisons.(comparison_name).(slot);
    full_S = project.models.(reference).model.S;
    
    fluxsum_cell = cellfun(@(x) get_model_fluxsum(full_S(:,x), samples(x,:),flux_summed_up), ...
                           rxn_idx,'UniformOutput', false);
    if flux_summed_up ~= "reactions"
        % this is only needed in case the fluxsum is computed over all the
        % metabolites in the given reactions, as an additional filtering
        % step, in case of "reactions" only one value is returend per set
        % for the get_model_fluxsum
        fluxsum_cell = cellfun(@(f) f(met_idx,:), ...
                               fluxsum_cell, 'UniformOutput', false);
    end

    

end
