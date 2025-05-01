function [sampling_ordered] = get_sampling_orig_order(m,s,num_rxns_orig)
                                sampling_values = zeros(num_rxns_orig,1000);
                                sampling_values(m.AA,:) = s;
                                sampling_ordered = sampling_values;
end