function fluxsum_cell = get_fluxsum(project,comparison_name,met_idx,rxn_idx,slot)
    arguments
        project
        comparison_name
        met_idx =[]
        rxn_idx =[]
        slot  (1,1) string {mustBeMember(slot,["ordered_fba", "ordered_samples"])} ="ordered_samples" 
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
    
    fluxsum_cell = cellfun(@(x) get_model_fluxsum(full_S(:,x), samples(x,:)), ...
                           rxn_idx,'UniformOutput', false);

    fluxsum_cell = cellfun(@(f) f(met_idx,:), ...
                           fluxsum_cell, 'UniformOutput', false);

    

end
