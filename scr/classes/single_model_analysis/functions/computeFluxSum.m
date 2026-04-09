function [fluxSum] = computeFluxSum(model, analysis, analysisType, fluxSumType, interestName)   
% model and analysis come from writeAnalysisReport
% analysisType = FBA/sampling
% fluxSumType = metabolite/pathway
% fluxSum = table
% interestName = name of the metabolite/pathway

%% Verifying the arguments
arguments
    model
    analysis
    analysisType {mustBeTextScalar}
    fluxSumType   {mustBeTextScalar}
    interestName
end

analysisType = string(analysisType);
fluxSumType   = string(fluxSumType);

mustBeMember(analysisType, ["FBA","sampling"])
mustBeMember(fluxSumType, ["metabolite","pathway"])

%% Get indices of rxns

if fluxSumType == "metabolite"
    pattern = ['^', interestName, '\['];
    idx = ~cellfun('isempty', regexp(model.mets, pattern, 'once'));
    matchedMets = model.mets(idx);
    [metsRxns, ~] = findRxnsFromMets(model, matchedMets);
    idxRxns = ismember(model.rxns, metsRxns); % indices of rxns associated with met of interest
    idxMets = ismember(model.mets, matchedMets); % indices in S of the mets of interest (several compartments)
else % fluxSumType = pathway
    % list of rxns associated with the pathway of interest
    pathwayRxns = findRxnsFromSubSystem(model, interestName);
    idxRxns = ismember(model.rxns, pathwayRxns);
end

%% Computing fluxSum

fluxSum = table();

if analysisType == "FBA"
    FBA = analysis.FBA.x;
    % Get S and FBA filtered with only the reactions of interest
    filteredS = model.S(:, idxRxns);
    filteredFBA = FBA(idxRxns);
    % Expanding
    expandedFilteredFBA = repmat(filteredFBA',size(filteredS,1),1);
    % Flux Sum
    fluxesFBA = filteredS.*expandedFilteredFBA;
    % Positive/Negative
    fluxesSumPos = full(sum((fluxesFBA>0).*fluxesFBA, 2)); % each element is multiplied by 1 if positive, otherwise it becomes 0 : matrix with positive flux only
    fluxesSumNeg = full(sum((fluxesFBA<0).*fluxesFBA, 2)); % opposite
    
    % Formatting the output table
    if fluxSumType == "metabolite"
        fluxSum.Metabolite = matchedMets;
        fluxSum.FluxSumPositiveFBA = fluxesSumPos(idxMets);
        fluxSum.FluxSumNegativeFBA = fluxesSumNeg(idxMets);
    else % fluxSumType = pathway
        fluxSum.Pathway = interestName;
        fluxSum.FluxSumPositiveFBA = sum(fluxesSumPos); % sum the rows = sum of the fluxes that go through all the reactions of the pathway
        fluxSum.FluxSumNegativeFBA = sum(fluxesSumNeg);
    end

else
    samples = analysis.sampling.samples;
    % Filter rxns
    filteredS = model.S(:, idxRxns);
    filteredSamples = samples(idxRxns);
    % Creation of the matrix which is going to store the flux sum of each sample
    allFluxesSumPos = zeros(size(filteredS,1),size(samples,2));
    allFluxesSumNeg = zeros(size(filteredS,1),size(samples,2));
    % Loop through the samples
    for i = 1:size(filteredSamples,2)
        sample = filteredSamples(:,i);
        % Expanding
        expandedSample = repmat(sample',size(filteredS,1),1);
        % Flux Sum
        fluxesSample = filteredS.*expandedSample;
        % Positive/Negative
        fluxesSumPos = full(sum((fluxesSample>0).*fluxesSample, 2)); % each element is multiplied by 1 if positive, otherwise it becomes 0 : matrix with positive flux only
        fluxesSumNeg = full(sum((fluxesSample<0).*fluxesSample, 2)); % opposite
        % Storing results
        allFluxesSumPos(:, i) = fluxesSumPos;
        allFluxesSumNeg(:, i) = fluxesSumNeg;
    end

    % Formatting the output table
    if fluxSumType == "metabolite"
        fluxSum.Metabolite = matchedMets;
        fluxSum.AverageFluxSumPositive = mean(allFluxesSumPos(idxMets), 2);
        fluxSum.MinFluxSumPositive = min(allFluxesSumPos(idxMets), [], 2);
        fluxSum.MaxFluxSumPositive = max(allFluxesSumPos(idxMets), [], 2);
        fluxSum.AverageFluxSumNegative = mean(allFluxesSumNeg(idxMets), 2);
        fluxSum.MinFluxSumNegative = min(allFluxesSumNeg(idxMets), [], 2);
        fluxSum.MaxFluxSumNegative = max(allFluxesSumNeg(idxMets), [], 2);
    else % fluxSumType = pathway
        fluxSum.Pathway = interestName;
        fluxSum.AverageFluxSumPositive = mean(sum(allFluxesSumPos, 1)); % sum the rows = sum of the fluxes that go through all the reactions of the pathway
        fluxSum.MinFluxSumPositive = min(sum(allFluxesSumPos, 1));
        fluxSum.MaxFluxSumPositive = max(sum(allFluxesSumPos, 1));
        fluxSum.AverageFluxSumNegative = mean(sum(allFluxesSumNeg, 1));
        fluxSum.MinFluxSumNegative = min(sum(allFluxesSumNeg, 1));
        fluxSum.MaxFluxSumNegative = max(sum(allFluxesSumNeg, 1));
    end
end
