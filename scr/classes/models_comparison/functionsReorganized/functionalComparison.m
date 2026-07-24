function plots = functionalComparison(project, comparisonName)
    % This function runs the functional model comparison. 
    % The models are compared on basis of the FBA & FVA results from the singleModelAnalysis. 
    % So the functional capacity the model has in context of the defined
    % objective function. 
    % Input:
    %   - project: struct with content defined in the README,
    %     singleModelAnalysis run on the object, and chooseActiveAnalysis
    %     needs to be set too
    %   - comparisonName: which of the comparisons should be visualized,
    %     comparison slot is created when running the
    %     modelStructuralComparison function
    % Output: #TODO
    %   - different figures
    %   - excel sheet - imports export fba, rxns wise FVA scores for each model, enrichment
    %     table with the scores -> needs to be done still 
    arguments
        project (1,1) struct
        comparisonName (1,1) string
    end

    modelList = project.comparisons.(comparisonName).modelNames;
    referenceModel = project.comparisons.(comparisonName).referenceModel;
    
    
    %%% ---------- Visualization: objective values per model
    fbaObjectiveValues = cell2mat(cellfun(@(x) project.models.(x).analysis.active.FBA.f(1,1), modelList, "UniformOutput", false));
    getExchangeRxnsIdx = find(findExcRxns(project.models.(referenceModel).model));    

    plots.objValue = figure('Color', 'w', 'Position', [20 20 700 300], 'Visible', 'off');
 
    bar(fbaObjectiveValues)
    title('Model comparison: flux of optimized reaction')
    ylabel('Reaction flux value for objective function [µMol/(gDW*h)]')
    xlabel('Model')
    xticklabels(modelList)
    set(gca, 'FontSize', 18)

    %%% ---------- Visualization: get FBA values under objective function
    %%%             for the different models - filtered for exchange rxns
    
    % Import
    plots.import = getFluxPlot(project, comparisonName, getExchangeRxnsIdx, ...
                                    'thresholdFlux', 'upper', 'FVA', false, 'reducedCost', false, 'visiblePlots', "off");
    % Export
    plots.export = getFluxPlot(project, comparisonName,getExchangeRxnsIdx,...
                  'thresholdFlux','lower','FVA',false,'reducedCost',false,'visiblePlots',"off");
    
    %%% ---------- Visualization: FVA Similarity between Models

    [fvaSim, fvaSimRxns, fvaSimPathways] = computeFvaSimilarity(project, comparisonName);

    plots.fvaSim.overall = plotClustergram(fvaSim,...
                             modelList,...
                             modelList,...
                             {'Similarity of FVA boundaries'},...
                             "FVA similarity",...
                             [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255);
    
    %%% ---------- Visualization: FVA Similarity per reaction histogramm

    plots.fvaSim.hist = FVASimValuesHist(fvaSimRxns, modelList);
    
    %%% ---------- Visualization: FVA Similarity per reaction - enrichment
    %%%            for low fva similarity scores per pathway in the model

    resEnrichment = getEnrichmentTable(project,modelList, fvaSimRxns, referenceModel,[]);
    % put the results of FDR and NES in one matrix each

    comparisons = fieldnames(resEnrichment);

    % All unique pathways
    allPathways = unique(vertcat(resEnrichment.(comparisons{1}).Subsystem));
    for k = 2:numel(comparisons)
        allPathways = unique([allPathways; resEnrichment.(comparisons{k}).Subsystem]);
    end

    % Preallocate tables
    NESTbl = array2table(nan(numel(allPathways), numel(comparisons)), ...
        'RowNames', allPathways, 'VariableNames', comparisons);
    FDRTbl = array2table(nan(numel(allPathways), numel(comparisons)), ...
        'RowNames', allPathways, 'VariableNames', comparisons);

    % Fill tables
    for c = 1:numel(comparisons)
        comp = comparisons{c};
        idx = ismember(allPathways, resEnrichment.(comp).Subsystem);
        [~, loc] = ismember(allPathways(idx), resEnrichment.(comp).Subsystem);
        NESTbl{idx, comp} = resEnrichment.(comp).NES(loc);
        FDRTbl{idx, comp} = resEnrichment.(comp).FDR(loc);
    end

    filterForSig = find((sum(FDRTbl{:,:} < 0.05,2) > 0) & (sum(NESTbl{:,:} > 0,2) >0));
    FDRTbl = FDRTbl(filterForSig,:);
    NESTbl = NESTbl(filterForSig,:);

    plots.fvaSim.enrich = dotplot(NESTbl,FDRTbl);
    
    
    %%% ---------- Visualization: Fluxsum based on the FBA values ? 

    replacementValue = "analysis.active.FBA.v"; % get the fba solution values
    project.comparisons.(comparisonName).orderedFba = getOrderedFeatureMatrix(project, modelList, "rxns", referenceModel, replacementValue);
    
    % compute Fluxsum 

    [idxPathways, names_pathways] = getDefaultSubsystems(project, referenceModel);                                 

    [fluxsumSets,plots.fba.heatmapRxnFluxsum] = visualizeFluxsum(project, comparisonName, [], idxPathways,...
                                                                     names_pathways,...
                                                                     "heatmap", true, true, "orderedFba", "reactions",...
                                                                      referenceModel, "off");
    
    
    plots.fba.heatmapRxnActivityFba = getNetworkActivity(project, comparisonName, idxPathways, names_pathways);


    [fluxsumSets,plots.fba.heatmapMetsFluxsum] = visualizeFluxsum(project, comparisonName, [], idxPathways,...
                                                                      names_pathways,...
                                                                      "heatmap", true, true, "orderedFba", "incoming",...
                                                                      referenceModel, "off");

    %%% -> show the top 20 most variant metabolites excluding known cofactors 
    % cofactorNames = ["atp", "adp", "amp", "nad", "nadh", "nadp", "nadph", ...
    %             "coa", "accoa", "fad", "fadh2", "pi", "pp_i"];

    %%% -> show the fluxsum in the defined pathways
    % get_essential_pathway_metabolites('Glycolysis',project,referenceModel)

    function fig = FVASimValuesHist(fvaSimRxns, modelList)
        % This function visualizes the FVA values in a histogramm per
        % comparison. These histogramms give us an indication of how similary
        % models are in their FVA min and max boundaries, and although the
        % heatmap summing up the FVA values also hast the same intuition, this
        % visualization also enables to see what similarity values are mostly
        % occuring. Is the difference in the overall similarity mainly due to a
        % few reactions having very low values, or do we see a lot of mean fva
        % similarity values per reaction ? 
        % #TODO improve function documentation!!
        fvaLower2x2 = getLowerTriangleBlock(fvaSimRxns);
        
        modelPairs = cell(numel(modelList));  % preallocate
        
        for i = 1:numel(modelList)
            for j = 1:numel(modelList)
                modelPairs{i,j} = {};
                if i ~= j
                    modelPairs{i,j} = {modelList{i}, modelList{j}};
                end
            end
        end
        modelPairs2x2 = getLowerTriangleBlock(modelPairs);
        
        [nRows, nCols] = size(fvaLower2x2);
        
        fig = figure('Color','w','Visible','off','Position', [100 100 2000*3 2000]);
        % Create tiled layout
        t = tiledlayout(fig,nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
        
        for i = 1:nRows
            for j = 1:nCols
                nexttile((i-1)*nCols + j)
        
                data = fvaLower2x2{i,j};
    
                if ~isempty(data)
                    data = data(data ~= 1);  % remove trivial values
                end
            
                if ~isempty(data)
                    histogram(data,100)
                    set(gca, 'FontSize', 18)
                else
                    axis off  % empty tile
                end
        
                % Label axes
                if ~isempty(data)
                    if i ~= 1
                        xlabel(modelPairs2x2{i,j}{1,2}, 'Interpreter','none')
                    end
                    if j == 1
                        ylabel(modelPairs2x2{i,j}{1,1}, 'Interpreter','none')
                    end
                end
                box on
            end
        end
        
        sgt = sgtitle('Histogram of FVA similarity values between models per rxns (<1)');
        sgt.FontSize = 20;   % set desired font size
    end

    function Results = getEnrichmentTable(project, modelList, fvaSimRxns, referenceModel, subSystems)
        % This function visualizes the enrichment results in a dotplot!!
        % #TODO: better documentation of the function1!!!
    
        [~,rxnMapping] = getOrderedFeatureMatrix(project, modelList, "rxns", referenceModel);
    
    
        if isempty(subSystems)
            subSystems = string(project.models.(referenceModel).model.subSystems); 
        end
        rxns = string(project.models.(referenceModel).model.rxns);        
    
        [uniqSubs, ~, idx] = unique(subSystems);
    
        subStruct = struct();
    
        for k = 1:numel(uniqSubs)
            subName = uniqSubs{k};
            fieldName = matlab.lang.makeValidName(subName);
        
            % Collect reactions belonging to this subsystem
            subStruct.(fieldName).name = subName;
            subStruct.(fieldName).rxns = rxns(idx == k);
        end
        n = length(modelList);
        [I, J] = ndgrid(1:n, 1:n);
        modelIndex = arrayfun(@(i,j) [i j], I, J, 'UniformOutput', false);
    
        fvaSim = getLowerTriangleBlock(fvaSimRxns);
        modelIndex = getLowerTriangleBlock(modelIndex);
     
        modelPairs = cell(numel(modelList));  % preallocate
        
        for i = 1:numel(modelList)
            for j = 1:numel(modelList)
                modelPairs{i,j} = {};
                if i ~= j
                    modelPairs{i,j} = {modelList{i}, modelList{j}};
                end
            end
        end
        modelPairs2x2 = getLowerTriangleBlock(modelPairs);
        
        
        Results    = struct();
        for k = 1:numel(fvaSim)
    
            x = fvaSim{k};
            y = strjoin(modelPairs2x2{k},'_');
    
            model1idx = modelIndex{k}(1);
            model2idx = modelIndex{k}(2);
    
            if isempty(x) || isempty(y)
                continue
            end
    
            rxnIdsInBothModels = find(sum(rxnMapping(:, [model1idx,model2idx]) ~= 0, 2) == 2);
            % filter for the rxn similarities that are in both models 
            rxnsInBothModels = rxns(rxnIdsInBothModels);
            
            Results.(string(y)) = pathwayEnrichment(subStruct, x(rxnIdsInBothModels), rxnsInBothModels);
    
        end
    
    
        function results = pathwayEnrichment(sets, metricMatrix, featureNames)
        % This function performs pathway enrichment on the fva similarity
        % values. In the context of metabolic modelling the enrichment in this
        % function is defined as the enrichment of low fva similarity values
        % (high dissimilarity) in a specific metabolic pathway. 
        % So practically this function does a ranked based hyptothesis testing.
        % Sorting of the rxns after their similarity value (ascending sorting,
        % since we want to know where the FVA boundaries are most different)
        % and then see which of the metabolic subsystems defined in the model
        % are enriched in this sorting!
        % # TODO -> better documentation of the function + check what happens
        % with the rnxs that have the same rank, there should be a group fo
        % rxns that have the same rxn FVA similarity -> does this effect the
        % enrichment -> since the sorting with the same value is then kind of
        % arbitraty !!!! 
    
    
        [metricSorted, sortIdx] = sort(metricMatrix, 'descend');
        rxnsSorted = featureNames(sortIdx);
        N = numel(featureNames);
        nPerm = 1000;   % permutations
        p = 1;          % weight exponent (0 = unweighted)
        weights = abs(metricSorted).^p;
    
    
        subNames = fieldnames(sets);
        nSets = numel(subNames);
        
        ES   = nan(nSets, 1);
        NES  = nan(nSets, 1);
        pval = nan(nSets, 1);
        setSize = nan(nSets, 1);
        
        for s = 1:nSets
        
            subField = subNames{s};
            subRxns  = sets.(subField).rxns;
            setSize(s) = numel(subRxns);
        
            % Skip very small subsystems
            if setSize(s) < 5
                continue
            end
        
            % Hits in ranked list
            hits = ismember(rxnsSorted, subRxns);
            Ns = sum(hits);
        
            % ----- observed enrichment score -----
            Phit  = weights .* hits;
            Phit  = Phit / sum(Phit);
            Pmiss = (~hits) / (N - Ns);
        
            runningSum = cumsum(Phit - Pmiss);
        
            [~, imax] = max(abs(runningSum));
            ES(s) = runningSum(imax);
        
            % ----- permutation test -----
            ESnull = zeros(nPerm, 1);
        
            for k = 1:nPerm
                permHits = hits(randperm(N));
        
                Phit_p  = weights .* permHits;
                Phit_p  = Phit_p / sum(Phit_p);
                Pmiss_p = (~permHits) / (N - Ns);
        
                rs_p = cumsum(Phit_p - Pmiss_p);
                ESnull(k) = max(abs(rs_p));
            end
        
            % empirical p-value
            pval(s) = mean(abs(ESnull) >= abs(ES(s)));
        
            % normalized enrichment score
            NES(s) = ES(s) / mean(abs(ESnull));
        end
    
        qval = mafdr(pval, 'BHFDR', true);
        results = table( ...
        subNames, setSize, ES, NES, pval, qval, ...
        'VariableNames', {'Subsystem', 'Size', 'ES', 'NES', 'pValue', 'FDR'} );
    
        results = sortrows(results, 'NES', 'descend');
    
    
        %results = results(results.FDR < 0.05, :);
    
    
        % Choose subsystem to plot
        % subField = 'HeparanSulfateDegradation'; 
        % subRxns  = sets.(subField).rxns;
        % 
        % % Compute hits
        % hits = ismember(rxnsSorted, subRxns);
        % Ns = sum(hits);
        % 
        % % Compute running sum for enrichment
        % Phit  = weights .* hits;
        % Phit  = Phit / sum(Phit);
        % Pmiss = (~hits) / (N - Ns);
        % 
        % runningSum = cumsum(Phit - Pmiss);
        % 
        % % Plot
        % figure('Color','w','Position',[200 200 600 400])
        % plot(runningSum, 'b', 'LineWidth', 2)
        % hold on
        % yline(0,'k--','LineWidth',1)
        % 
        % % Mark hits as vertical lines along x-axis
        % hitIdx = find(hits);
        % for i = 1:numel(hitIdx)
        %     x = hitIdx(i);
        %     line([x x], [min(runningSum) 0], 'Color',[0.7 0.7 0.7],'LineStyle','-')
        % end
        % 
        % xlabel('Reactions ranked by FVA similarity')
        % ylabel('Running enrichment score')
        % title(['Enrichment of ', sets.(subField).name])
        % grid on
    
    
    end
 end   
function fig = dotplot(NESTbl, FDRTbl)
    % #TODO better documentation of the function!!
    
    % --- Sort pathways by overall NES magnitude ---
    [~, sortedIdx] = sort(sum(abs(NESTbl{:,:}), 2), 'descend');
    NESTbl = NESTbl(sortedIdx, :);
    FDRTbl = FDRTbl(sortedIdx, :);
    
    % --- Handle zeros in FDR and transform ---
    lowValues = 1e-10;
    FDRTbl{:,:}(FDRTbl{:,:} == 0) = lowValues;
    FDRTbl{:,:} = -log10(FDRTbl{:,:});
    
    % --- Extract labels ---
    pathways = regexprep(string(NESTbl.Properties.RowNames), "_", " ");
    comparisons = string(NESTbl.Properties.VariableNames);
    
    % --- Extract numeric matrices ---
    NES = NESTbl{:,:};
    FDR = FDRTbl{:,:};
    
    nP = numel(pathways);
    nC = numel(comparisons);
    
    % --- Create grid for scatter ---
    [X, Y] = meshgrid(1:nC, 1:nP);
    x = X(:);
    y = Y(:);
    
    nesVals = NES(:);
    fdrVals = FDR(:);
    
    % --- Prepare FDR for coloring ---
    cVals = fdrVals;              
    cVals(cVals < -log10(0.05)) = NaN;  % values >0.05 will be grey

    % --- Dot size proportional to |NES| with enhanced visual difference ---
    scatterMin = 10;    % smallest dot area (points^2)
    scatterMax = 500;  % largest dot area (points^2)

    nes = nesVals(~isnan(cVals));
    
    absNESNorm = (abs(nes) - min(abs(nes))) / (max(abs(nes)) - min(abs(nes))); % normalize 0-1
    dotSize = scatterMin + (absNESNorm.^0.5) * (scatterMax - scatterMin);  % power 0.5 emphasizes large values
    
    
    % --- Create figure ---
    fig = figure('Color', 'w', 'Position', [100 100 1000 1000], "Visible","off");
    hold on
    
    % Scatter plot for significant FDR (≤ 0.05)
    scatter(x(~isnan(cVals)), y(~isnan(cVals)), dotSize, cVals(~isnan(cVals)), 'filled')
    hold on
    
    % --- Size legend for |NES| with min, percentiles, max ---
    sizeVals = [min(abs(nes)), ...
                 prctile(abs(nes), 25), ...
                 prctile(abs(nes), 50), ...
                 prctile(abs(nes), 75), ...
                 max(abs(nes))];
    
    sizeScaled = [min(dotSize), ...
                 prctile(dotSize, 25), ...
                 prctile(dotSize, 50), ...
                 prctile(dotSize, 75), ...
                 max(dotSize)];
   
    % Custom "legend" inside axes
    legendSizes = sizeScaled;      % sizeData
    legendLabels = string(round(sizeVals, 2));
    legendX = max(x) + 1;  % x position outside plot
    legendY = y(1:length(legendSizes)) ;        % y positions
    
    for i = 1:length(legendSizes)
        scatter(legendX, legendY(i), legendSizes(i), 'k', 'filled')
        text(legendX+0.2, legendY(i), legendLabels{i}, 'FontSize', 12, 'VerticalAlignment','middle')
    end
    % --- Axes formatting ---
    xticks(1:nC)
    xticklabels(regexprep(comparisons, "_", " vs "))
    yticks(1:nP)
    yticklabels(pathways)
    
    xlabel('Model comparison')
    ylabel('Pathway')
    
    xlim([0.5, nC + 1.5])
    ylim([0.5, nP + 0.5])
    set(gca, 'YDir', 'reverse', 'FontSize', 18)
    title("Pathway enrichment (dot size = |NES|, color = FDR)")
    
    % --- Colorbar ---
    nColors = 256;
    cmap = [linspace(1,0,nColors)' linspace(0,0,nColors)' linspace(0, 1 , nColors)']; 
    colormap(cmap)        % red (low) -> blue (high)
    clim([-log10(0.05) - log10(lowValues)])       
    cb = colorbar;
    cb.Label.String = '-log10(FDR)';
    cb.FontSize = 14;
    
end

function fvaFower = getLowerTriangleBlock(fvaSimRxns)
        % This function gets the lower part of a similarity matrix. 
        % Written because we are looking at pairwise distances/similarities in
        % the figures, but to not have repetitive plots it is useful to only
        % look at one part of the matrix (eiter above or below the diagonal)
        % since the diagonal (so the comparison of the model with itself) will
        % always be 1 or zero (depending wheterh we talk about distance or
        % similarity) and therefore only one of the triangels in the matrix
        % cell is interesting. 
        % This function sets all but the lower triangle to 0 so that there are
        % no repetitive plots!
        % fvaSimRxns: n x n cell array of comparisons
        % Returns a compact lower-triangle block cell array
    
        n = size(fvaSimRxns, 1);
        if n ~= size(fvaSimRxns, 2)
            error('Input must be a square cell array');
        end
    
        % Count how many rows/columns to keep: remove first row if needed
        % General: keep lower-triangle below the first row
        rowsToKeep = 2:n;       % always skip the first row
        colsToKeep = 1:(n-1);   % always skip the last column
    
        % Take the lower-triangle block
        fvaFower = fvaSimRxns(rowsToKeep, colsToKeep);
end

end