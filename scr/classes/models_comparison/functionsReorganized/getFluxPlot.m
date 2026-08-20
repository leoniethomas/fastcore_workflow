function fig = getFluxPlot(project, comparisonName, idxToVis, options)
    % This function visualizes the values of FVA,FBA, and reduced Cost for
    % the choosen rxns_ids between different models.
    % As an input a singlemodelanalysis needs to be
    % given, the name of the comparison for which this plot is to be
    % displayed, and the indices of the rxns to be displayed. In the
    % options you can choose which values you want to be displayed, by
    % default only the FBA values are displayed in the barplot. But in
    % the options it can additionally be specified that the FVA boundaries
    % (as a grey box) and the reduced cost per reaction (as color of the
    % FBA dot) can be displayed. 
    % Input: 
    % - project:                project object which is the output of
    %                           single_model_analysis script
    % - comparison_name:        name of the comparison as a string
    % - idxToVis:             positions of the rxns to be displayed in the
    %                           choosen reference model
    % - options:                - FVA = true (default false) to display the FVA boundaries
    %                             around the FBA solution as a grey box.                        
    %                           - thresholdFlux= wether to apply an upper
    %                             lower ,no upper or lower or to include 
    %                             also the fba =0 reactions(all) to the selected
    %                             reaction fba values
    % 
    % Output:                   Display of figure                
    % 

    arguments
        project
        comparisonName (1,1) string
        idxToVis
        options.FVA  (1,1) logical = false
        %options.reducedCost (1,1) logical = false
        options.thresholdFlux (1,1) string {mustBeMember(options.thresholdFlux, ["lower", "upper", "none", "all"])} = "none" 
        options.titlePlots = ""
        options.visiblePlots = "on"
    end


    modelList = project.comparisons.(comparisonName).modelNames;
    referenceModel = project.comparisons.(comparisonName).referenceModel;
    
    replacementValue = "analysis.active.FBA.v"; % get the fba solution values
    orderedFbaMatrix = getOrderedFeatureMatrix(project, modelList, referenceModel, "rxns", replacementValue);
    replacementValue = "mappedDiscretizedRxns"; % get the fba solution values
    orderedMappingRxnMatrix = getOrderedFeatureMatrix(project, modelList, referenceModel, "rxns", replacementValue);
    
    if options.thresholdFlux == "upper"
        getExchangeRxnsIdx = intersect(find(sum(orderedFbaMatrix, 2) ~=0 & sum(orderedFbaMatrix < 0, 2) ~= 0), ...
                                          idxToVis);  
        titleWord = "Negative flux reactions";
    elseif options.thresholdFlux == "lower"
        getExchangeRxnsIdx = intersect(find(sum(orderedFbaMatrix,2) ~=0 & sum(orderedFbaMatrix > 0,2) ~= 0), ...
                                                      idxToVis); 
        titleWord = "Positive flux reactions";
    elseif options.thresholdFlux == "none"
        getExchangeRxnsIdx = intersect(find(sum(orderedFbaMatrix, 2) ~=0), ...
                                                      idxToVis);
        titleWord = options.titlePlots;
    elseif options.thresholdFlux == "all"
        getExchangeRxnsIdx = idxToVis;
        titleWord = options.titlePlots;
    else
        error("Wrong value for chosen threshold. Possible values: 'lower', 'upper', 'none', 'all'.")
    end

    orderedFbaMatrixEx = orderedFbaMatrix(getExchangeRxnsIdx, :);
    orderedMappingRxnMatrixEx = orderedMappingRxnMatrix(getExchangeRxnsIdx, :);
    rxnNames = project.models.(referenceModel).model.rxns(getExchangeRxnsIdx);

    if options.FVA
        replacementValue = "analysis.active.FVA.minMaxFluxes.maxFlux"; % get the fba solution values
        orderedFvaMaxMatrix = getOrderedFeatureMatrix(project, modelList, referenceModel, "rxns", replacementValue);
        replacementValue = "analysis.active.FVA.minMaxFluxes.minFlux"; % get the fba solution values
        orderedFvaMinMatrix = getOrderedFeatureMatrix(project, modelList, referenceModel, "rxns", replacementValue);
        orderedFvaMaxMatrixEx = orderedFvaMaxMatrix(getExchangeRxnsIdx, :);
        orderedFvaMinMatrixEx = orderedFvaMinMatrix(getExchangeRxnsIdx, :);
        
        if options.thresholdFlux == "all"
            idxFVAVar = ~(all(orderedFvaMaxMatrixEx == 0, 2) & all(orderedFvaMaxMatrixEx == 0, 2));
            orderedFvaMaxMatrixEx = orderedFvaMaxMatrixEx(idxFVAVar, :);
            orderedFvaMinMatrixEx = orderedFvaMinMatrixEx(idxFVAVar, :);
            rxnNames = rxnNames(idxFVAVar);
            orderedFbaMatrixEx = orderedFbaMatrixEx(idxFVAVar, :);
            orderedMappingRxnMatrixEx = orderedMappingRxnMatrixEx(idxFVAVar, :);
        end
        
        % if options.reducedCost
        % 
        %     models = rmfield(project.models,setdiff(fieldnames(project.models), modelList));
        % 
        %     if all(structfun(@(x) isfield(x.analysis.active.FBA, "w"), models)) && all(structfun(@(x) ismember('one', x.analysis.active.parameters.Value(x.analysis.active.parameters.Analysis == "FBA")), models))
        % 
        %         for m = modelList'
        %             nRxns = length(project.models.(m).model.rxns);
        %             w = project.models.(m).analysis.active.FBA.w;  % 7350×1
        %             project.models.(m).analysis.active.FBA.reducedCost = w(2*nRxns+1 : 3*nRxns);    % net flux variables ← use this one
        %         end
        % 
        %         replacementValue = "analysis.active.FBA.reducedCost"; % get the fba solution values
        %         orderedReducedCostMatrix = getOrderedFeatureMatrix(project, modelList, referenceModel, "rxns", replacementValue);
        %         orderedReducedCostMatrixEx = orderedReducedCostMatrix(getExchangeRxnsIdx, :);
        %         if options.thresholdFlux == "all"
        %             orderedReducedCostMatrixEx = orderedReducedCostMatrixEx(idxFVAVar, :);
        %         end
        %     else
        %         error("The reduced costs could not be found in the .w slot of the FBA analysis. In case you did not use the 'one' minNorm parameter for the FBA, shadowprices might be stored elsewere.")
        %     end
        % end
        
    end
 
    rxnFormulas = printRxnFormula(project.models.(referenceModel).model, 'rxnAbbrList', rxnNames, 'printFlag', false);

    % medium composition 
    % do the models have the same  medium ? 
    % if so I can add it as column to the plot, otherwise it is not added
    models = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    mediaForModels = structfun(@(x) x.settings.medium, models);
    ref = fieldnames(models);
    ref = models.(ref{1,1}).settings.medium;
    mediumIsEqualBetweenModels = all(arrayfun(@(x) isequaln(x, ref), mediaForModels));

    if mediumIsEqualBetweenModels
        if isfield(ref,'mediumComposition') %in case there is no medium defined
            rxnNames = string(rxnNames);
            % checking which column of the medium composition table entails
            % the rxn names -> just checking for the column with most
            % overlap
            isText = varfun(@(x) iscell(x) || isstring(x) || ischar(x), ...
                            ref.mediumComposition, ...
                            "OutputFormat", "uniform");

            textTable = ref.mediumComposition(:, isText);
            nMatches = varfun(@(x) sum(ismember(string(x), rxnNames)), textTable);
            [maxMatches, bestColumn] = max(nMatches{:,:});
            bestColumnName = ref.mediumComposition.Properties.VariableNames{bestColumn};

            mediumConstrained = ismember(rxnNames, ref.mediumComposition.(bestColumn)); %| ...
                             % ismember(rxns_names, ref.manual_set_boundaries.unwanted_export) | ...
                             % ismember(rxns_names, ref.manual_set_boundaries.unwanted_import);
        end

        orderedLb = getOrderedFeatureMatrix(project, referenceModel, ...
                                             referenceModel, "rxns", "model.lb");
        orderedUb = getOrderedFeatureMatrix(project, referenceModel, ...
                                             referenceModel, "rxns", "model.ub");
        orderedUb = orderedUb(getExchangeRxnsIdx, :);
        orderedLb = orderedLb(getExchangeRxnsIdx, :);
        
        if options.FVA && options.thresholdFlux == "all"
            orderedLb = orderedLb(idxFVAVar, :);
            orderedUb = orderedUb(idxFVAVar, :);
        end

        % get rxn gene rules to add to the table
        symbolGprRules = string(cellfun(@(rxnName)getRxnSymbolRule(project.models.(referenceModel), ...
                                                       rxnName), string(rxnNames), 'UniformOutput', false));
       
        T = table(rxnFormulas, mediumConstrained,...
            join(string(orderedMappingRxnMatrixEx), "|", 2), symbolGprRules, ...
                  'VariableNames', ["Reaction Formula", "Medium Constrained", ...
                                    join(string(project.comparisons.(comparisonName).modelNames), "_"), ...
                                    "Symbol GPR Rules"], 'RowNames', rxnNames);
        T = T(flip(string(T.Properties.RowNames)), :);
        T.lb = flip(orderedLb);
        T.ub = flip(orderedUb);
    end
    
    % if options.shadowPrice
    %     replacementValue = "analysis.active.FBA.basis.dual"; % get the fba solution values
    %     % shadow prices are measured for every metabolite therefore mapped according to the mets field
    %     ordered_shadowPrices_matrix = getOrderedFeatureMatrix(project, modelList, referenceModel, "mets", replacementValue);
    % end

    %%% Parameters for the figure 

    %%% =========================
    % Threshold & masks
    % =========================
    mu = mean(orderedFbaMatrixEx, 2);
    
    x  = abs(mu);
    
    % 1D k-means to minimize within-group distances

    assert(length(x) >1,...
           'Only one reaction value was returned after filtering for the active reactions in the FBA solution. In order to visualize this plot, you need to formulate more rxns in rxnID or to apply other threholds (try "all" threshold).')
    
    [idx, C] = kmeans(x, 2, 'Replicates', 10);
    
    % Identify low-flux (dense) cluster
    [~, lowCluster] = min(C);
    
    isLow = idx == lowCluster;
    isHigh = ~isLow;
    
    % Derive lowRange for documentation
    % lowRange = [-2 2];
    lowRange = [-max(x(isLow)), max(x(isLow))];
    
    isLow  = mu > lowRange(1) & mu < lowRange(2);
    isHigh = ~isLow;
    
    dataLow  = orderedFbaMatrixEx(isLow, :);
    dataHigh = orderedFbaMatrixEx(isHigh, :);
    
    rxnNamesLow  = rxnNames(isLow);
    rxnNamesHigh = rxnNames(isHigh);
    
    %%% =========================
    % Adaptive relative heights
    % =========================
    alpha = 0.75;
    nLow  = numel(rxnNamesLow);
    nHigh = numel(rxnNamesHigh);
    
    w = [nLow nHigh].^alpha;
    usableHeight = 1 - 0.08 - 0.10 - 0.03;
    
    heightLow  = usableHeight * w(1) / sum(w);
    heightHigh = usableHeight * w(2) / sum(w);
    
    bottomLow  = 0.10;
    bottomHigh = bottomLow + heightLow + 0.03;
    %%%
    
    fig = uifigure('Name', titleWord + " with Table", ...
                       'Position', [100 100 1000 450], 'Visible', options.visiblePlots);
        
    plotWidth = 0.52;
    
    %%% =========================
    % Helper for axes formatting
    % =========================
    fmtAxes = @(ax) set(ax, ...
        'FontSize', 16, ...
        'PositionConstraint', 'outerposition', ...
        'GridColor', [.8 .8 .8], ...
        'GridAlpha', .5);
    
    if ~options.FVA %&& ~options.reducedCost % specify which of the fields need to be true and false!!!
        % in the case that only the FBA solution should be visualized we
        % use a grouped horizontal barplot to do so

        %%% =========================
        % TOP AXIS — high values
        % =========================
        if heightHigh ~= 0
            axTop = uiaxes(fig, 'Units', 'normalized', ...
                'Position', [0.02 bottomHigh plotWidth heightHigh]);
            fmtAxes(axTop)

            if size(dataHigh, 1) == 1
                dataHigh = [dataHigh; NaN(1, size(dataHigh, 2))];
                dataHighOneBar = 1;
            else
                dataHighOneBar = 0;
            end
            
            barh(axTop, dataHigh, 'grouped');
            grid(axTop, 'on')
            
            %axTop.YTick = [];
            %axTop.YColor = 'none';
            
            yticks(axTop, 1:numel(rxnNamesHigh))
            yticklabels(axTop, strrep(rxnNamesHigh, "_", "\_"))
            title(axTop, titleWord)
            if heightLow == 0
                xlabel(axTop, 'Flux value [µMol/(gDW*h)]')
            end
            if dataHighOneBar
                ylim(axTop, [0.5, size(dataHigh, 1) - 0.4 ])
            end
        end

        %%% =========================
        % BOTTOM AXIS — low values
        % =========================
        if heightLow ~= 0
            axBottom = uiaxes(fig, 'Units', 'normalized', ...
                'Position', [0.02 bottomLow plotWidth heightLow]);
            fmtAxes(axBottom)
            
            if size(dataLow, 1) == 1
                dataLow = [dataLow; NaN(1, size(dataLow, 2))];
                dataLowOneBar = 1; 
            else
                dataLowOneBar = 0;
            end
            barh(axBottom, dataLow, 'grouped');
            grid(axBottom, 'on')
            
            yticks(axBottom, 1:numel(rxnNamesLow))
            yticklabels(axBottom, strrep(rxnNamesLow, "_", "\_"))
            if heightHigh == 0
                title(axBottom, titleWord)
            end
            xlabel(axBottom, 'Flux value [µMol/(gDW*h)]')
            if dataLowOneBar
                ylim(axBottom, [0.5, size(dataLow, 1) - 0.4 ])
            end
        end
        
        %%% =========================
        % Reverse direction if needed
        % =========================
        if options.thresholdFlux == "upper"
            set([axTop axBottom], 'XDir', 'reverse')
        end
        
        %%% =========================
        % Legend
        if nHigh < nLow
            legend(axBottom, modelList, 'Location', 'best');
        else
            legend(axTop, modelList, 'Location', 'best');
        end

    else

        Q1_low  = orderedFvaMinMatrixEx(isLow, :);
        MED_low = orderedFbaMatrixEx(isLow, :);
        Q3_low = orderedFvaMaxMatrixEx(isLow, :);
        
        Q1_high  = orderedFvaMinMatrixEx(isHigh, :);
        MED_high = orderedFbaMatrixEx(isHigh, :);
        Q3_high = orderedFvaMaxMatrixEx(isHigh, :);
            
        % =========================
        % Check if reduced cost coloring is enabled
        % =========================
        %addColorLegend = options.reducedCost;
        
        % if addColorLegend
        %     cmap = colormap('cool');  % N colors
        %     N = size(cmap, 1);
        % 
        %     % Split reduced cost into high and low
        %     reducedCostLow  = orderedReducedCostMatrixEx(isLow, :);
        %     reducedCostHigh = orderedReducedCostMatrixEx(isHigh, :);
        % 
        %     % Helper function to scale matrix to colormap indices
        %     scaleToCmap = @(mat, valMin, valMax) round((mat - valMin) / (valMax - valMin) * (N-1)) + 1;
        % 
        %     % Get global min/max from full matrix
        %     valMin = min(orderedReducedCostMatrixEx(:));
        %     valMax = max(orderedReducedCostMatrixEx(:));
        % 
        %     if valMin == valMax
        %         % Entire matrix has one unique value — all same color
        %         scaledIdxLow = ones(size(reducedCostLow));
        %         scaledIdxHigh = ones(size(reducedCostHigh));
        %     else
        %         % Check each subset individually
        %         if length(unique(reducedCostLow(:))) == 1
        %             scaledIdxLow = ones(size(reducedCostLow));
        %         else
        %             scaledIdxLow = scaleToCmap(reducedCostLow, valMin, valMax);
        %         end
        % 
        %         if length(unique(reducedCostHigh(:))) == 1
        %             scaledIdxHigh = ones(size(reducedCostHigh));
        %         else
        %             scaledIdxHigh = scaleToCmap(reducedCostHigh, valMin, valMax);
        %         end
        %     end
        % else
        %     addColorLegend = 0;
        % end

       %%
        
        % --- TOP AXES: high flux ---
        %%% =========================
        % Plot high reactions
        %%% =========================
        if nHigh > 0
            axTop = uiaxes(fig, 'Units', 'normalized', 'Position', [0.02 bottomHigh plotWidth heightHigh]);
            hold(axTop, 'on'); 
            fmtAxes(axTop);

            % --- Grid & font ---
            grid(axTop, 'on')
            axTop.GridColor = [0.8 0.8 0.8];
            axTop.GridAlpha = 0.5;
            axTop.FontSize = 16;
        
            nGroups = nHigh;
            nPerGroup = size(Q1_high, 2);
            boxHeight = 0.2; 
            groupSep = 1;
            greyLevels = linspace(0.9, 0.6, nPerGroup);
            greyRGBs = [greyLevels', greyLevels', greyLevels'];
        
            for i = 1:nGroups
                yBase = i*groupSep;
                for j = 1:nPerGroup
                    offset = (j-(nPerGroup+1)/2)*(boxHeight*1.3);
                    rectangle(axTop, 'Position', [Q1_high(i,j), yBase+offset-boxHeight/2, Q3_high(i,j)-Q1_high(i,j), boxHeight],...
                              'FaceColor', greyRGBs(j, :), 'EdgeColor', 'none');
                    % Median dot
                    % if addColorLegend
                    %     plot(axTop, MED_high(i, j), yBase+offset, 'o', ...
                    %     'MarkerSize', 5, ...
                    %     'MarkerFaceColor', cmap(scaledIdxHigh(i, j), :), ...
                    %     'MarkerEdgeColor', cmap(scaledIdxHigh(i, j), :));
                    % else
                        plot(axTop, MED_high(i, j), yBase+offset, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
                    %end
                end
            end
        
            yticks(axTop, 1:nGroups*groupSep)
            yticklabels(axTop, strrep(rxnNamesHigh, "_", "\_"))
            title(axTop, titleWord)
        end
        
        %%% =========================
        % Plot low reactions
        %%% =========================
        if nLow > 0
            
            axBottom = uiaxes(fig, 'Units', 'normalized', 'Position', [0.02 bottomLow plotWidth heightLow]);
            hold(axBottom, 'on'); 
            fmtAxes(axBottom);
            
            % --- Grid & font ---
            grid(axBottom, 'on')
            axBottom.GridColor = [0.8 0.8 0.8];
            axBottom.GridAlpha = 0.5;
            axBottom.FontSize = 16;
        
            nGroups = nLow;
            nPerGroup = size(Q1_low, 2);
            boxHeight = 0.2; 
            groupSep = 1;
            greyLevels = linspace(0.9, 0.6, nPerGroup);
            greyRGBs = [greyLevels', greyLevels', greyLevels'];
            
            for i = 1:nGroups
                yBase = i*groupSep;
                for j = 1:nPerGroup
                    offset = (j-(nPerGroup+1)/2)*(boxHeight*1.3);
                    rectangle(axBottom, 'Position', [Q1_low(i,j), yBase+offset-boxHeight/2, Q3_low(i,j)-Q1_low(i,j), boxHeight], ...
                              'FaceColor', greyRGBs(j,:), 'EdgeColor', 'none');
                    % Median dot
                    % if addColorLegend
                    %     plot(axBottom, MED_low(i,j), yBase+offset, 'o', ...
                    %     'MarkerSize', 5, ...
                    %     'MarkerFaceColor', cmap(scaledIdxLow(i, j), :), ...
                    %     'MarkerEdgeColor', cmap(scaledIdxLow(i, j), :));
                    % else
                        plot(axBottom, MED_low(i, j), yBase+offset, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
                %     end
                end
            end

            yticks(axBottom, 1:nGroups*groupSep)
            yticklabels(axBottom, strrep(rxnNamesLow, "_", "\_"))
            xlabel(axBottom, 'Flux value [µMol/(gDW*h)]')
        end

        %%% =========================
        % Set axes limits based on MED only
        %%% =========================
        if exist('axTop', 'var')
            % X limits
            xMinHigh = min(MED_high(:));
            xMaxHigh = max(MED_high(:));
            % Add small padding
            xPad = (xMaxHigh - xMinHigh) * 0.05;
            % xlim(axTop, [xMinHigh-xPad, xMaxHigh+xPad]);
        
            % Y limits
            yMinHigh = 0.5;  % first y
            yMaxHigh = nHigh * groupSep + 0.5;
            ylim(axTop, [yMinHigh, yMaxHigh]);
        end
        
        if exist('axBottom','var')
            % X limits
            xMinLow = min(MED_low(:));
            xMaxLow = max(MED_low(:));
            xPad = (xMaxLow - xMinLow)*0.05;
            if xMinLow-xPad == xMaxLow-xPad
                xlim(axBottom, [-1, 1]);
            else
                xlim(axBottom, [xMinLow-xPad, xMaxLow+xPad]);
            end
        
            % Y limits
            yMinLow = 0.5;
            yMaxLow = nLow * groupSep + 0.5;
            ylim(axBottom, [yMinLow, yMaxLow]);
        end
        
        % --- Add colorbar if needed
        % if addColorLegend
        %     % Choose axTop if it exists, else axBottom
        %     if exist('axTop', 'var')
        %         cbAx = axTop;
        %     else
        %         cbAx = axBottom;
        %     end
        % 
        %     colormap(cbAx, cmap);          % Set colormap for chosen axis
        %     % --- Set caxis safely
        %     if valMin == valMax
        %         caxis(cbAx, [valMin-0.5, valMax+0.5]);  % small padding if single value
        %     else
        %         caxis(cbAx, [valMin, valMax]);
        %     end
        % 
        %     cb = colorbar(cbAx);           % Attach colorbar
        %     cb.Label.String = ...
        %        'Reduced Cost';
        %     cb.FontSize = 14;
        % end
        
        % --- Grey patches legend (for models)
        % Choose axTop if it exists, else axBottom
        if exist('axTop', 'var')
            lgdAx = axTop;
        else
            lgdAx = axBottom;
        end
        
        hGrey = gobjects(nPerGroup, 1);
        for j = 1:nPerGroup
            hGrey(j) = patch(lgdAx, NaN, NaN, greyRGBs(j,:), 'EdgeColor', 'none');
        end
        
        % Add legend
        lgd = legend(lgdAx, hGrey, modelList, 'Location', 'northeastoutside');
        % lgd = legend(lgdAx, hGrey, modelList, 'Location','northeast');
        lgd.FontSize = 16;
        lgd.Box = 'off';
        lgd.Color = 'none';
        lgd.Title.String = "Models";

    end

    % --- Reverse direction if needed ---
    if options.thresholdFlux == "upper"
        set([axTop axBottom], 'XDir', 'reverse')
    end
    
    %%% =========================
    % Table
    % =========================
    % resort the table according to the order in the plots
    if exist('T', 'var')
        T = T([flip(string(rxnNamesHigh)); flip(string(rxnNamesLow))], :);
        tbl = uitable(fig, ...
                        'Data', T, ...
                        'ColumnName', T.Properties.VariableNames, ...
                        'Units', 'normalized', ...
                        'Position', [plotWidth+0.05 0.10 0.40 0.85], ... % width increased from 0.30 → 0.40
                        'FontSize', 16, ...
                        'ColumnWidth', 'auto');
         tbl.FontSize = 16; 
    end

end




