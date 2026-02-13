function fluxsum = get_model_fluxsum(S_matrix,samples,flux_summed_up)
    arguments
        S_matrix 
        samples 
        flux_summed_up {mustBeMember(flux_summed_up, ["incoming","outgoing"])} ="incoming" 
    end
        
        % get rid of zero entries
        non_zero_rxns = find(~all(samples == 0,2));
        % subset S and sampling matrix accordingly 
        S_matrix = S_matrix(:,non_zero_rxns);
        samples = samples(non_zero_rxns,:);
        
        fluxsum=zeros(size(S_matrix,1),size(samples,2));
        for counter=1:size(samples,2)
            v=samples(:,counter); % one sample
            temp=repmat(v',size(S_matrix,1),1); %
            fluxes=S_matrix.*temp;
            
            fluxSumP=full(sum((fluxes>0).*fluxes,2));
            fluxSumN=full(sum((fluxes<0).*fluxes,2));

            if abs(fluxSumP) ~= abs(fluxSumN)
                error("Your fluxes seem to be missassigned. There was some mix up, resulting in the incoming Fluxsum not being equal to the outgoing fluxsum. ")
            end
            if flux_summed_up == "outgoing"
                fluxsum(:,counter)=fluxSumN;
            else
                fluxsum(:,counter)=fluxSumP;
            end
        end
end