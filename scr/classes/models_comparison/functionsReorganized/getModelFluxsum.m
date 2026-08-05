function fluxsum = getModelFluxsum(S_matrix, samples, fluxSummedUp)
    % This functions computes the fluxsum over all solutions given in the
    % samples matrix. Fluxes can be summed up either for the metabolites or
    % for a set of reactions given. 
    % So the output fluxsum is therefore either a measure of usage of a set
    % of reactions or a measure of usage of one or more metabolites (for
    % the metabolite either the incoming or the outgoing rxns are summed
    % up)
    % INPUT: 
    %   - S_matrix:         stochiometric matrix of the model from which the
    %                       solution(s) come from given in samples matrix. Only needed for the
    %                       metabolite fluxsum, for the reaction fluxsum it can be left empty
    %                       (default: []).
    %   - samples:          matrix with as many reaction flux values as there are
    %                       reactions in the stochiometric matrix
    %   - fluxSummedUp:   gives the fluxes over which to sum up, either
    %                       computation of the fluxsum for the single
    %                       metabolites ("incoming" or "outgoing") or
    %                       summing up the flux values over the defined
    %                       number of reactions: "reactions". 
    %                       "Incoming" and "outgoing" just defines whether
    %                       the positive or the negative fluxes connected
    %                       to a metabolite will be summed up (abs sum).
    % OUTPUT:
    %   - fluxsum: fluxsum for the defined reactions, either fluxsum over
    %   the metabolites or over the rxns defined, the output is the size of
    %   the metabolites in the model in case "incoming" or "outgoing" is
    %   choosen, and only one value for the "reactions" option.
    
    arguments
        S_matrix 
        samples 
        fluxSummedUp {mustBeMember(fluxSummedUp, ["incoming","outgoing","reactions"])} ="incoming" 
    end
        
    % check that stochiometric matrix is the same size as samples
    % matrix
    if size(S_matrix, 2) ~= size(samples, 1)
        error("The stochiometric matrix and the samples do not have the same count of rxns. It seems the stochiometric matrix is not derived from the model for which the solutions were drawn.")
    end

    % get rid of zero entries
    nonZeroRxns = find(~all(samples == 0, 2));
    % subset S and sampling matrix accordingly 
    S_matrix = S_matrix(:, nonZeroRxns);
    samples = samples(nonZeroRxns, :);

    if fluxSummedUp == "reactions"
        % fluxsum over the reaction fluxes, for that we do not need
        % stochiometry
        fluxsum = sum(abs(samples), 1); % better to sum over absolute fluxes ?   
    else
        % fluxsum for each metabolite, either all outgoing or all incoming
        
        fluxsum = zeros(size(S_matrix, 1), size(samples, 2));
        for counter = 1:size(samples, 2)
            v = samples(:, counter); % one sample
            temp = repmat(v', size(S_matrix, 1), 1); %
            fluxes = S_matrix.*temp;
            
            fluxSumP = full(sum((fluxes > 0).*fluxes, 2));
            fluxSumN = full(sum((fluxes < 0).*fluxes, 2));

            if abs(fluxSumP) ~= abs(fluxSumN)
                error("Your fluxes seem to be missassigned. There was some mix up, resulting in the incoming Fluxsum not being equal to the outgoing fluxsum.")
            end
            if fluxSummedUp == "outgoing"
                fluxsum(:, counter) = fluxSumN;
            else
                fluxsum(:, counter) = fluxSumP;
            end
        end
    end
end