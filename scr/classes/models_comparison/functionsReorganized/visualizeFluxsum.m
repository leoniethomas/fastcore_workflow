function [fluxsumSets, fig] = visualizeFluxsum(project, comparisonName, metIdx, rxnIdx, rxnSetLabels, plotType, excludeCoenzymes, ignoreCompartment, slot, fluxSummedUp, modelNameFluxBoundaries, plotVisible)
    % This function visualizes the fluxsum either in a heatmap or in a
    % violinPlot between different models to enable easy comparison of
    % model results in sampling. 
    % Before running this script a modelsComparison must be run. Then the
    % name of the comparison performed (given as output of the
    % modelsComparison function) can be given as an input together with the
    % project object entailing the comparison + all the needed models. 
    %
    % The function produces as many figures as you define sets in the
    % rxnIdx argument. This is meant to enable you to visualize different
    % sets/subsystems of your choosing to be visualized in one Figure. 
    % 
    % The function allows visualization of the fluxsum on different levels.
    % If you choose the violin plotType -> the fluxsum for every metabolite is
    % computed and then visualized in a separate violinPlot for every
    % reaction allowing for easy comparison of the fluxsum of over
    % different models. 
    % If you choose the heatmap then the mean of the fluxsum overall
    % metabolites in the defined rxn_set will be computed (excluding the
    % coenzymes if specified) and therefore one value per model per rxn_id
    % set is can be visualized in a heatmap to give a more overall sense of
    % the activity of for example different subsystems in the model. 
    % 
    % INPUT:
    %       - project: project structure after running modelsComparison
    %         function 
    %       - comparisonName: name of the comparison you want to visualize
    %         (the second output of the modelsComparison function) 1x1
    %         string
    %       - metIdx: array 1xn, every element in the array is a double
    %         1xn with the positions of the mets to be used in the fluxsum
    %         comptutation, the positions need to be the position of that
    %         metabolite in the choosen reference model
    %       - rxnIdx: same as metIdx
    %       - rxnSetLabels: names of the defined sets, can be choosen
    %         freely
    %       - plotType: to visualize a heatmap (fluxsum overall given in
    %         one rnx_id set) or violin (fluxsum for all metabolites taking
    %         part in the specified rxns)
    %       - excludeCoenzymes: option to exclude all the metabolite that
    %         have a connectivity of over 135 within the reference network
    %         (those are most likely coenzymes -> like h2o)
    %       - ignoreCompartment: allows to ignore the compartment and
    %         build an overall fluxsum per metabolite (only relevant for
    %         violin plots)
    %       - fluxSummedUp: gives the fluxes over which to sum up, either
    %                       computation of the fluxsum for the single
    %                       metabolites ("incoming" or "outgoing") or
    %                       summing up the flux values over the defined
    %                       number of reactions: "reactions". 
    %                       "Incoming" and "outgoing" just defines whether
    %                       the positive or the negative fluxes connected
    %                       to a metabolite will be summed up (abs sum).
    % Output: 
    %       - fluxsumSets: fluxsum overall metabolites per sample as a
    %       array with matrices stored within
    arguments
        project struct
        comparisonName (1,1) string
        metIdx % (1,:) cell {mustBeColumnVector} = {}
        rxnIdx (1,:) cell {mustBeColumnVector} = {}
        rxnSetLabels (1,:) string = ""
        plotType {mustBeMember(plotType, ["violin", "heatmap", "heatmapSample", "heatmapSampleAllFeatures"])} = ["violin"] 
        excludeCoenzymes (1,1) logical = false
        ignoreCompartment (1,1) logical = true
        slot  (1,1) string {mustBeMember(slot, ["orderedFba", "orderedSamples","orderedllSamples"])} = "orderedSamples" 
        fluxSummedUp {mustBeMember(fluxSummedUp, ["incoming", "outgoing", "reactions"])} = "incoming" 
        modelNameFluxBoundaries (1,1) string = "consistentMediumConstrainedModel"
        plotVisible ="on"
    end
    
    % sets the order of the rxns in order to be able to bring the samples
    % from all the models in the same order
    reference = project.comparisons.(comparisonName).referenceModel;

    % by default if no input is given the fluxsum for all metabolites overa
    % ll reaction is computed
    if isempty(rxnIdx)
        rxnIdx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    if isempty(metIdx)
        metIdx = find(ones(length(project.models.(reference).model.mets),1));
    end

    if isstring(metIdx)
        metIdx = find(cellfun(@(m) any(strcmp(regexp(m, '[^\[]+', 'match', 'once'), cellstr(string(metIdx)))), project.models.(reference).model.mets));
    end
    
    if plotType == "heatmapSample" & slot == "orderedFba"
        % if you choose to plot the fba solutions, the heatmapSample
        % option does not exist, cause there is only one solution per model
        plotType = "heatmap";
    end
    
    % in order to prevent coenzymes from inflating the fluxsum of specific
    % pathways all metabolites which have a connectivity higher than 135 in
    % the choosen reference model will be excluded from the computation of
    % the fluxsum and not be visualized in the violinPlot
    if excludeCoenzymes
        S = project.models.(reference).model.S;
        rxnCountPerMet = sum(S ~= 0, 2);  % sum over columns = number of reactions each metabolite participates in
        idxCoenzymes = find(rxnCountPerMet >135);

        % histogram(rxn_count, 1000)
        % xlim([0,50])
        % xlabel('Number of reactions per metabolite')
        % ylabel('Count of metabolites')
        
        % define coenzymes by connectivity in the network, highly connected
        % metabolites are probably coenzymes 
        % idx_coenzymes = find(rxn_count_per_met > quantile(rxn_count_per_met,0.994));

        metIdx = setdiff(metIdx, idxCoenzymes);
    end
    
    % depending on the choosen plot type different functions are executed
    if plotType == "violin"
            [fluxsumSets, fig] = getViolinPlots(project, comparisonName, metIdx, rxnIdx, rxnSetLabels, ignoreCompartment, modelNameFluxBoundaries, "incoming",slot, plotVisible);
    else
            [fluxsumSets, fig] = getComparisonHeatmap(project, comparisonName, metIdx, rxnIdx, rxnSetLabels, plotType, slot, fluxSummedUp, plotVisible);
    end

    % when metIdx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites
end

function [fluxsumSets, fig] = getComparisonHeatmap(project, comparisonName, metIdx, rxnIdx, rxnSetLabels, type, slot, fluxSummedUp, plotVisible)
    % This function visualizes the mean fluxsum over the specified rxn_id
    % sets. By getting the fluxsum for metabolites that are participating
    % in the specified sets + are part of the metIdx and then computing
    % the mean over each set. Using this function a comparison over models
    % and different defined subsystems can be visualized. 
    % INPUT:
    %       - project: project structure after running modelsComparison
    %         function 
    %       - comparisonName: name of the comparison you want to visualize
    %         (the second output of the modelsComparison function) 1x1
    %         string
    %       - metIdx: with the positions of the mets to be used in the fluxsum
    %         comptutation, the positions need to be the position of that
    %         metabolite in the choosen reference model
    %       - rxnIdx: sets of rxns to be visualized in the heatmap
    %       - rxnSetLabels: names of the defined sets, can be choosen
    %         freely
    %       - type: specifies whether a fluxsum per rxn set should be
    %         visualized per samples or the average over all samples from one
    %         model: values : "heatmapSample", "heatmapSampleAllFeatures" or "heatmap" default is "heatmap"
    %       - fluxSummedUp:   gives the fluxes over which to sum up, either
    %                       computation of the fluxsum for the single
    %                       metabolites ("incoming" or "outgoing") or
    %                       summing up the flux values over the defined
    %                       number of reactions: "reactions". 
    %                       "Incoming" and "outgoing" just defines whether
    %                       the positive or the negative fluxes connected
    %                       to a metabolite will be summed up (abs sum).
    % Output: 
    %       - fluxsumSets: fluxsum overall metabolites per sample as a
    %         array with matrices stored within
    arguments
        project struct
        comparisonName (1,1) string
        metIdx  
        rxnIdx
        rxnSetLabels (1,:) string 
        type {mustBeMember(type, ["heatmap", "heatmapSample", "heatmapSampleAllFeatures"])} = ["heatmap"] 
        slot  (1,1) string {mustBeMember(slot, ["orderedFba", "orderedSamples", "orderedllSamples"])} = "orderedSamples" 
        fluxSummedUp {mustBeMember(fluxSummedUp, ["incoming", "outgoing", "reactions"])} = "incoming" 
        plotVisible = "on"
    end
    
    reference = project.comparisons.(comparisonName).referenceModel;
    
    if isempty(rxnIdx)
        rxnIdx = find(ones(length(project.models.(reference).model.rxns), 1));
    end
    
    if isempty(metIdx)
        metIdx = find(ones(length(project.models.(reference).model.mets), 1));
    end
    
    if slot == "orderedFba" & type == "reactions" || type == "heatmapSampleAllFeatures"
        % this case is an exception
        % in this case we just visualize the flux values - no fluxsum
        % calculation needed 
        fbaValues = project.comparisons.(comparisonName).(slot);
        fluxsumSets = cellfun(@(x) fbaValues(x,:), ...
                           rxnIdx, 'UniformOutput', false);
    else
        % for all the other cases we compute the fluxsum 
        fluxsumSets = getFluxsum(project, comparisonName, metIdx, rxnIdx, slot, fluxSummedUp);
    end

    if slot ~= "orderedFba"
        samplesCat = cellstr(project.comparisons.(comparisonName).sampleModelLabels);
    else
        samplesCat = cellstr(project.comparisons.(comparisonName).modelNames);
    end

    if fluxSummedUp ~= "reactions"
        figureTitle = "FluxSum (per metabolite) ";
    else
        figureTitle = "FluxSum (over defined reaction set) ";
    end
    if slot ~= "orderedFba"
        figureTitle = figureTitle + " for sampling solutions ";
    else
        figureTitle = figureTitle + " for FBA solution ";
    end

    heatmapData = zeros(length(fluxsumSets), length(unique(samplesCat)));
    heatmapDataAllSamples = zeros(length(fluxsumSets), length(samplesCat));

    % when metIdx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites
    heatmapDataAllSamplesAllFeatures = {};
    for subsystem = 1:numel(fluxsumSets)
        data = fluxsumSets{subsystem};
        titleFig = rxnSetLabels(subsystem);
        
        metNames = project.models.(reference).model.mets(metIdx);
    
        metNames = metNames(find(any(data ~= 0, 2)));
        data = data(find(any(data ~= 0, 2)), :);

        %
        samplesCat = categorical(samplesCat); 
        groups = unique(samplesCat, 'stable');
        nGroups = numel(groups);
        [nMet, nSamples] = size(data);
        dataGrouped = cell(1,nGroups);
    
        for g = 1:nGroups
            idx = samplesCat == groups(g);   % logical index for this group
            dataGrouped{g} = data(:, idx);   % all metabolites, only this group
        end

        heatmapData(subsystem, :) = cellfun(@(x) mean(x(:)), dataGrouped);

        heatmapData(isnan(heatmapData)) = 0;

        heatmapDataAllSamples(subsystem,:) = cell2mat(cellfun(@(x) mean(x,1), dataGrouped, 'UniformOutput', false));
        heatmapDataAllSamplesAllFeatures{end +1} = cell2mat(dataGrouped);

    end

    if type == "heatmap"
        % z-scaling for the heatmap in order to make the differences between
        % samples for one pathway more visible
        scaledData = zscore(heatmapData')';
        
        fig = figure('Color', 'w', 'Position', [100 100 800 800], 'Visible', plotVisible);
        imagesc(scaledData)
        
        cmap = getColorPallette();
        h = colorbar;  
        caxis([-max([max(scaledData(:)), abs(min(scaledData(:)))]) max([max(scaledData(:)), abs(min(scaledData(:)))])])   % colors scaled from -2 (min) to 2 (max)

        ylabel(h, 'Scaled average fluxsum average value over rxn set', 'FontSize', 18)        % Set title/label of colorbar
        axis equal tight % Make cells square and remove extra space
    
        title(figureTitle + "average overall solutions") % grayscale
        % Set x-axis and y-axis labels
        set(gca, 'XTick', 1:length(unique(samplesCat)), 'XTickLabel', unique(samplesCat,'stable'), ...
             'YTick', 1:length(rxnSetLabels), 'YTickLabel', rxnSetLabels)
        xtickangle(45)
        ax = gca;
        ax.FontSize = 18;  
        xlabel('Model', 'FontSize', 18)       
        ylabel('Reaction set', 'FontSize', 18)    
       
        [nRows, nCols] = size(heatmapData);
    
        % Loop over every cell and place the absolute number from pathway_counts
        for i = 1:nRows
            for j = 1:nCols
                % You want absolute numbers, not relative counts, so use pathway_counts
                % (or multiply relative_counts by referenceModel if needed)
                value =  heatmapData(i, j); % +1 because first column is referenceModel
                % Place text at the center of the tile
                text(j, i, num2str(adaptiveFormat(value)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Color', 'k', ...          % black text
                    'FontSize', 18)
            end
        end
        
        hold off
    
    elseif type == "heatmapSample"

        fig = figure('Color', 'w', 'Position', [100 100 800 800], 'Visible', plotVisible);
        scaledData = zscore(heatmapDataAllSamples')';
    
        imagesc(scaledData)

        cmap = getColorPallette();
       
        h = colorbar;  
        caxis([-max([max(scaledData(:)), abs(min(scaledData(:)))]) max([max(scaledData(:)), abs(min(scaledData(:)))])])   % colors scaled from -2 (min) to 2 (max)

        ylabel(h, 'Scaled average fluxsum value per sample', 'FontSize', 18) % Set title/label of colorbar
        title(figureTitle + " values for each solution") % grayscale
        % Set x-axis and y-axis labels
        sampleCount = histcounts(samplesCat); % or your existing "a"

        edges = [0 cumsum(sampleCount)];
        xtickposition = edges(1:end-1) + sampleCount/2;
    
        set(gca, 'XTick', xtickposition, 'XTickLabel', unique(samplesCat,'stable'), ...
             'YTick', 1:length(rxnSetLabels), 'YTickLabel', rxnSetLabels)
        xtickangle(45)
        ax = gca;
        ax.FontSize = 18;  
        xlabel('Model', 'FontSize', 18)       
        ylabel('Reaction set', 'FontSize', 18)
    
    else

    % z-scaling for the heatmap in order to make the differences between
        % samples for one pathway more visible!
        scaledData = zscore(cell2mat(heatmapDataAllSamplesAllFeatures')')';
        
        fig = figure('Color', 'w', 'Position', [100 100 800 800], 'Visible', plotVisible);
        imagesc(scaledData)
        
        cmap = getColorPallette();
        h = colorbar;  
        minAxis = quantile(scaledData(:), 0.001);
        maxAxis = quantile(scaledData(:), 0.999);
        axisLimit = max([abs(minAxis), abs(maxAxis)]);
        caxis([-axisLimit, axisLimit]) 

        ylabel(h, 'Scaled average fluxsum average value in reaction set', 'FontSize', 18)        % Set title/label of colorbar
        
        title(figureTitle + "value for each solution and feature") % grayscale
        % Set x-axis and y-axis labels
        [sampleCount, ~] = hist(samplesCat);
        xtickposition = ((1:length(unique(samplesCat))) .* (sampleCount)) - sampleCount/2;
    
        countMetsPerSet = cell2mat(arrayfun(@(x)size(x{:}, 1), heatmapDataAllSamplesAllFeatures, 'UniformOutput', false))';
        ytickposition = cumsum(countMetsPerSet) - round(countMetsPerSet/2);
    
        set(gca, 'XTick', xtickposition, 'XTickLabel', unique(samplesCat, 'stable'), ...
             'YTick', ytickposition, 'YTickLabel', rxnSetLabels)
        xtickangle(45)
        ax = gca;
        ax.FontSize = 18;  
        xlabel('Model', 'FontSize', 18)       
        ylabel('Reaction set', 'FontSize', 18)    
    
    end

end

function [fluxsumSets, figs] = getViolinPlots(project, comparisonName, metIdx, rxnIdx, rxnSetLabels, ignoreCompartment, modelNameFluxBoundaries, fluxSummedUp, slot, plotVisible)
    % This function visualizes the fluxsum distribution per metabolite and model 
    % over the specified rxn_id sets into a violin plot.
    % By getting the fluxsum for metabolites that are participating
    % in the specified sets + are part of the metIdx.
    % Using this function a comparison over models for the fluxsum of
    % different metabolites can be visualized.  
    % INPUT:
    %       - project: project structure after running modelsComparison
    %         function 
    %       - comparisonName: name of the comparison you want to visualize
    %         (the second output of the modelsComparison function) 1x1
    %         string
    %       - metIdx: with the positions of the mets to be used in the fluxsum
    %         comptutation, the positions need to be the position of that
    %         metabolite in the choosen reference model
    %       - rxnIdx: rxn sets, for each set a separate figure is created
    %         with all the metabolite fluxsum values for all the metabolite
    %         participating in the specified rxnIdx
    %       - rxnSetLabels: names of the defined sets, can be choosen
    %         freely
    %       - ignoreCompartment: allows to compute the fluxsum for a
    %         metabolite ignoring the cellular location
    %       - modelNameFluxBoundaries: model that saves the
    %         constraints set when constructing the models
    %   - fluxSummedUp:   gives the fluxes over which to sum up, either
    %                       computation of the fluxsum for the single
    %                       metabolites ("incoming" or "outgoing") or
    %                       summing up the flux values over the defined
    %                       number of reactions: "reactions". 
    %                       "Incoming" and "outgoing" just defines whether
    %                       the positive or the negative fluxes connected
    %                       to a metabolite will be summed up (abs sum).
    % Output: 
    %       - fluxsumSets: fluxsum overall metabolites per sample as a
    %         array with matrices stored within
    arguments
        project struct
        comparisonName (1,1) string
        metIdx
        rxnIdx
        rxnSetLabels (1,:)
        ignoreCompartment (1,1) logical = true
        modelNameFluxBoundaries (1,1) string = "consistentMediumConstrainedModel"
        fluxSummedUp {mustBeMember(fluxSummedUp, ["incoming", "outgoing", "reactions"])} = "incoming" 
        slot (1,1) string {mustBeMember(slot, ["orderedSamples","orderedllSamples"])} = "orderedSamples" 
        plotVisible = "on"
    end

    reference = project.comparisons.(comparisonName).referenceModel;
    referenceModel = project.models.(reference).model;
    modelNames = project.comparisons.(comparisonName).modelNames;
    models = rmfield(project.models, setdiff(fieldnames(project.models), modelNames));

    if isempty(rxnIdx)
        rxnIdx = find(ones(referenceModel.rxns), 1);
    end
    
    if isempty(metIdx)
        metIdx = find(ones(length(referenceModel.mets), 1));
    end

    fluxsumSets = getFluxsum(project, comparisonName, metIdx, rxnIdx,slot, fluxSummedUp);

    % for every specified set in rxnIdx create one figure with all the
    % fluxsums of the metabolites participating in the rxns + being part of
    % metIdx
    figs = struct();
    for subsystem = 1:numel(fluxsumSets)
        data = fluxsumSets{subsystem};
        titleFig = rxnSetLabels(subsystem);
        plotName = replace(titleFig, ["_", "-", "/"], "");

        metNames = referenceModel.mets(metIdx); % filter for mets in metIdx 
        
        zeroFluxSumMetabolites = find(any(data ~= 0, 2)); % filter for metabolites which have zero fluxsum overall samples, overall models
        data = data(zeroFluxSumMetabolites, :);
        metNames = metNames(zeroFluxSumMetabolites);
        samplesCat = cellstr(project.comparisons.(comparisonName).sampleModelLabels);
        % remove the compartment specification to compute the totall
        % fluxsum of a metabolite, and add them up overall compartments
        if ignoreCompartment
            metNamesWithoutCompartment = regexprep(metNames, "\[.\]$", "");

            % Get unique specifications and grouping indices
            [uniqueMets, ~, idx] = unique(metNamesWithoutCompartment);
            
            % Prepare output
            numSpecs  = numel(uniqueMets);
            numCols   = size(data, 2); 
            aggMatrix = zeros(numSpecs, numCols);
            
            % Sum rows by group
            for i = 1:numSpecs
                aggData(i, :) = sum(data(idx == i, :), 1);
            end

            data = aggData;
            metNames = uniqueMets;

        end
        
        samplesCat = categorical(samplesCat);  
        groups = unique(samplesCat, 'stable');       
        nGroups = numel(groups);
        [nMet, nSamples] = size(data);
        dataGrouped = cell(1,nGroups);
        
        % put the samples from different models in different matrices
        for g = 1:nGroups
            idx = samplesCat == groups(g);   
            dataGrouped{g} = data(:, idx);   
        end
        
        % put each metabolite in a different vector
        for m = 1:nMet
            for g = 1:nGroups
                dataForPlot{m, g} = dataGrouped{g}(m, :);  % 1 × nSamples_in_group
            end
        end

        % now we have a separate vector for every metabolite for every model
        % now we loop over all metabolites to visualize on violin per model
        
        maxPlotsPerFig = 12;
        nFigs = ceil(nMet/maxPlotsPerFig);
        
        figStruct = struct();

        for f = 1:nFigs
            
            fig = figure('Color', 'w', 'Position', [100 100 800 800], 'Visible', plotVisible);
            
            t = tiledlayout(3, 4);
            title(t,"Fluxsum for metabolites in : " + titleFig, ...
                'FontSize', 20, 'FontWeight','bold', 'Interpreter','none');
            
            % Store using dynamic field name
            fieldName = sprintf('fig%d', f);
            figStruct.(fieldName) = fig;
            
            startIdx = (f-1)*maxPlotsPerFig + 1;
            endIdx = min(f*maxPlotsPerFig, nMet);
            
            for m = startIdx:endIdx
                
                ax = nexttile(t);
                hold(ax, 'on')
                
                dat = dataForPlot(m, :);
                % in order to be able to put it all in one matrix we need
                % to make the sample count the same length, so we add NaNs
                
                maxLen = max(cellfun(@numel, dat));
                dat = cell2mat(cellfun(@(x) [x, nan(1, maxLen-numel(x))], dat, 'UniformOutput', false));
                dat = reshape(dat, maxLen, []);
                columnsToKeep = find(~all(dat == 0 | isnan(dat)));
                
                evalc('violinplot(dat(:, columnsToKeep), groups(columnsToKeep), "ShowData", false);');
                
                ylabel(ax, 'Flux Value', 'FontSize', 18);
                title(ax, metNames{m}, 'Interpreter','none', 'FontSize', 14);
                
                ax.FontSize = 18;
            end
        end
        
        if plotVisible == "on"
            % get Table with rxn formula, discretization status, concentration
            % in the medium etc
            warning('off', 'all')
            T_display = getRxnOverviewTable(project, comparisonName, models, referenceModel, metNames, modelNameFluxBoundaries, modelNames, reference);
            warning('on', 'all')
    
            T_display = addvars(T_display, string(T_display.Properties.RowNames), ...
                                'Before', 1, 'NewVariableNames', "Reaction");
            
            % Create UI figure
            fig = uifigure('Name', "Metabolite Fluxsum in subSystem: " + titleFig, ...
                           'Position', [100 100 1200 600], 'Visible', plotVisible);
            
            % Create table filling the entire figure
            tbl = uitable(fig, ...
                'Data', T_display{:,:}, ...
                'ColumnName', T_display.Properties.VariableNames, ...
                'Units', 'normalized', ...
                'Position', [0 0 1 1], ...
                'FontSize', 18, ...
                'ColumnWidth', 'auto');
        end
        
        plotName = replace(plotName, ["_", "-", "/", " "], "");
        figs.("violinFluxSum" + plotName) = figStruct;

    end

end

function T_display = getRxnOverviewTable(project, comparisonName, models, referenceModel, metNames, modelNameFluxBoundaries, modelNames, reference)
        
        ref = fieldnames(models);
        ref = models.(ref{1,1}).settings.medium;

        if all(~ismember(metNames, referenceModel.mets))
            % this happens when we ignore the compartment, cause then the
            % metabolites are without their compartment specification [c] 
            % in this case we search for biggest overlap
            metNamesIdx = find(cellfun(@(m) any(strcmp(regexp(m, '[^\[]+', 'match', 'once'), cellstr(string(metNames)))), referenceModel.mets));
            metNames = string(referenceModel.mets(metNamesIdx));
        end
        
        % get rxnIdx
        [rxnNames] = findRxnsFromMets(referenceModel, string(metNames));
        rxnIds = find(matches(referenceModel.rxns, rxnNames));
         
        % get rxnIdx formulas + filter out rxns that have zero in all
        % samples
        samples = project.comparisons.(comparisonName).orderedSamples;
        zeroRxns = find(sum(samples == 0, 2) == size(samples, 2));
        rxnIds = setdiff(rxnIds, zeroRxns);
        rxnNames = referenceModel.rxns(rxnIds);
        
        % get lower and upper bound for every rxns 
        orderedLb = getOrderedFeatureMatrix(project, modelNameFluxBoundaries, reference, "rxns", "model.lb");
        orderedUb = getOrderedFeatureMatrix(project, modelNameFluxBoundaries, reference, "rxns", "model.ub");
        orderedMappingRxnMatrix = getOrderedFeatureMatrix(project, modelNames, reference, "rxns", "mappedDiscretizedRxns");
    
        % check with samplings are not zero 
        if isfield(ref,'mediumComposition') %in case there is no medium defined
            rxnsNames = string(rxnNames);
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

        orderedUb = orderedUb(rxnIds, :);
        orderedLb = orderedLb(rxnIds, :);
        orderedMappingRxnMatrix = orderedMappingRxnMatrix(rxnIds, :);
        rxnAbbr = referenceModel.rxns(rxnIds);
        rxnFormulas = string(printRxnFormula(referenceModel, rxnAbbr, false));
  
        % get rxn gene rules to add to the table
        symbolGPRrules = string(cellfun(@(rxnName)getRxnSymbolRule(project.models.(reference), ...
            rxnName), string(rxnNames), 'UniformOutput', false));

        T = table(rxnFormulas, mediumConstrained, ...
            join(string(orderedMappingRxnMatrix), "|", 2),symbolGPRrules,...
                  'VariableNames', ["Reaction Formula","Medium Constrained", ...
                                    join(string(project.comparisons.(comparisonName).modelNames), "_"), ...
                                    "Symbol GPR Rules"], 'RowNames', rxnNames);
        T = T(flip(string(T.Properties.RowNames)), :);
        T.lb = flip(orderedLb);
        T.ub = flip(orderedUb);

        % Convert row names into a column
        T_display = T;

end

function mustBeColumnVector(c)
    % control function that checks that the rxnid, metid sets used as an input
    % in the visualize_fluxsum function are in the correct format 
    % array { 1xn double , 1xn double, 1xn double}
    for k = 1:numel(c)
        if ~isvector(c{k}) || size(c{k},2) ~= 1
            error('Each cell element must be an n×1 column vector.')
        end
    end
end

function str = adaptiveFormat(val)
% Format numbers adaptively:
% small values (|x| < 0.01) → scientific notation with 1 decimal places
% large values → fixed 1 decimal places

    if val == 0
        str = '0';
    elseif abs(val) < 0.01
        % Format as X.XX × 10^N
        exp = floor(log10(abs(val)));
        mantissa = val / 10^exp;
        str = sprintf('%.1f×10^{%d}', mantissa, exp);
    else
        str = sprintf('%.1f', val);
    end
end



