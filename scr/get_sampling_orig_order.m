function [sampling_ordered] = get_sampling_orig_order(m,s,rxns_orig)
                                [~,mapping_rxns_in_orig_idx] = ismember(m.rxns,rxns_orig);
                                sampling_values = zeros(length(rxns_orig),size(s,2));
                                sampling_values(mapping_rxns_in_orig_idx,:) = s;
                                id_biomass_ordered = find(matches(rxns_orig,"biomass_reaction"));
                                id_biomass = find(matches(m.rxns,"biomass_reaction"));
                                
%                                 min(sampling_values(id_biomass_ordered,:))
%                                 max(sampling_values(id_biomass_ordered,:))
%                                 
%                                 min(s(id_biomass,:))
%                                 max(s(id_biomass,:))
                                sampling_ordered = sampling_values;
end