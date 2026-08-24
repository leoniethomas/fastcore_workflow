function [fluxSets, figs] = visualizeFlux(project, comparisonName, rxnIdx, rxnSetLabels, slot, threshold, plotVisible,kld)
    arguments
        project
        comparisonName
        rxnIdx =[]
        rxnSetLabels = []
        slot (1,1) string {mustBeMember(slot, ["orderedSamples","orderedllSamples"])} = "orderedSamples" 
        threshold {mustBeMember(threshold, ["all", "positive", "negative"])} = ["all"] 
        plotVisible = "on"
        kld = 0
    end
    
    %     plot_type  {mustBeMember(plot_type, ["violin","heatmap"])} =["violin"] 
    %     exclude_coenzymes (1,1) logical = true
    %     ignore_compartment (1,1) logical = true
    % end
    
    reference = project.comparisons.(comparisonName).referenceModel;
    
    if isempty(rxnIdx)
        rxnIdx = {find(ones(length(project.models.(reference).model.rxns), 1))};
    end
    
    [fluxSets, figs] = getViolinPlotsFlux(project, comparisonName, rxnIdx, rxnSetLabels, threshold, slot, plotVisible,kld);

end

function [fluxSets, figs] = getViolinPlotsFlux(project, comparisonName, rxnsIdx, rxnSetLabels, threshold, slot,plotVisible,kld)
    
    arguments
        project
        comparisonName
        rxnsIdx =[]
        rxnSetLabels = []
        threshold {mustBeMember(threshold, ["all", "positive", "negative"])} = ["positive"] 
        slot (1,1) string {mustBeMember(slot, ["orderedSamples","orderedllSamples"])} = "orderedSamples" 
        plotVisible = "on"
        kld = 0
    end

    reference = project.comparisons.(comparisonName).referenceModel;
    referenceModel = project.models.(reference).model;
    modelNames = project.comparisons.(comparisonName).modelNames;
    models = rmfield(project.models, setdiff(fieldnames(project.models), modelNames));
    
    if isempty(rxnsIdx)
        rxnsIdx = {find(ones(length(project.models.(reference).model.rxns), 1))};
    end

    fluxSets = getFlux(project, comparisonName, rxnsIdx,slot);

    if kld
        KLD_sig = project.comparisons.(comparisonName).samplingComparison.kld.intraModelKld(cell2mat(rxnsIdx),:); % get KLD values
        KLD_sigLabels = project.comparisons.(comparisonName).samplingComparison.kld.intraModelKldLabels;
    end
    % when metIdx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites

    figs = struct();

    for subsystem = 1:numel(fluxSets)
        data = fluxSets{subsystem};
        if kld
            dataKLD = KLD_sig;
        end
        titleFig = string(rxnSetLabels(subsystem));
        plotName = replace(titleFig, ["_", "-", "/"], "");

        rxnsNames = project.models.(reference).model.rxns(rxnsIdx{subsystem});

        zeroFluxSumRxns = find(any(data ~= 0, 2));
        data = data(zeroFluxSumRxns, :);
        if  kld
            dataKLD = dataKLD(zeroFluxSumRxns, :);
        end
        rxnsNames = rxnsNames(zeroFluxSumRxns);

        if threshold ~= "all"
            Q = zeros(128, 4);
            for i = 1:size(data, 1)
                Q(i, :) = quantile(data(i, :), [0.20, 0.40, 0.60, 0.80]);
            end
            % compute quantiles -> check that at least 3  out of the 4
            % quantiles are either positive or negative

            if threshold == "positive"
                idx = find(sum(Q > 0, 2) > 2);
            elseif threshold == "negative"
                idx = find(sum(Q < 0, 2) > 2);
            end

            data = data(idx, :);
            if kld
                dataKLD = dataKLD(idx, :);
            end
            rxnsNames = rxnsNames(idx);
        end

        samplesCat = project.comparisons.(comparisonName).samplingComparison.sampleModelLabels;
        groups = unique(samplesCat, 'stable');       
        nGroups = numel(groups);
        [nMet, nSamples] = size(data);
        dataGrouped = cell(1, nGroups);
        
        for g = 1:nGroups
            idx = samplesCat == groups(g);   % logical index for this group
            dataGrouped{g} = data(:, idx);   % all metabolites, only this group
        end
    
        for m = 1:nMet
            for g = 1:nGroups
                dataForPlot{m,g} = dataGrouped{g}(m,:);  % 1 × nSamples_in_group
            end
        end
        
        maxPlotsPerFig = 12;
        nFigs = ceil(nMet/maxPlotsPerFig);
        
        figStruct = struct();

        for f = 1:nFigs
            
            fig = figure('Color', 'w', 'Position', [100 100 800 800], 'Visible', plotVisible);
            
            t = tiledlayout(3, 4);
            title(t,"Flux for reactions in : " + titleFig, ...
                  'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none');
            
            % Store using dynamic field name
            fieldName = sprintf('fig%d', f);
            figStruct.(fieldName) = fig;
            
            startIdx = (f-1)*maxPlotsPerFig + 1;
            endIdx = min(f*maxPlotsPerFig, nMet);
            
            for m = startIdx:endIdx
                
                ax = nexttile(t);
                hold(ax, 'on')
                
                dat = dataForPlot(m, :);
                if kld
                    datKLD = dataKLD(m, :);
                end
                % in order to be able to put it all in one matrix we need
                % to make the sample count the same length, so we add NaNs
                
                maxLen = max(cellfun(@numel, dat));
                dat = cell2mat(cellfun(@(x) [x, nan(1, maxLen-numel(x))], dat, 'UniformOutput', false));
                dat = reshape(dat, maxLen, []);
                columnsToKeep = find(~all(dat == 0 | isnan(dat)));
                evalc('violinplot(dat(:, columnsToKeep), groups(columnsToKeep), "ShowData", false);');

                hold on
                if kld
                    datKLD = reorderKLDComparisons(groups(columnsToKeep), datKLD, KLD_sigLabels);
                    yposition_text = max(max(dat(:, columnsToKeep)));
                    addSignificanceBars(datKLD , KLD_sigLabels,groups(columnsToKeep),dat(:, columnsToKeep))
                end
                
                % labels in plot
                groups = ["Control", "StageI", "StageII"];

                KLD_sigLabelsInPlot = strings(nchoosek(numel(groups), 2), 1);
                
                k = 1;
                for i = 1:numel(groups)-1
                    for j = i+1:numel(groups)
                        KLD_sigLabelsInPlot(k) = groups(i) + "_vs_" + groups(j);
                        k = k + 1;
                    end
                end

                %reorder datKLD accordingly 
                
                ylabel(ax, 'Flux Value', 'FontSize', 18);
                title(ax, rxnsNames{m}, 'Interpreter', 'none', 'FontSize', 14);
                
                ax.FontSize = 18;

            end
        end
        
        if plotVisible == "on"

        %table with rxnforumlas etc 

        ref = fieldnames(models);
        ref = models.(ref{1, 1}).settings.medium;
        %media_for_models = structfun(@(x) x.settings.medium, models);
        %medium_is_equal_between_models = all(arrayfun(@(x) isequaln(x, ref), media_for_models));

        rxnIds = find(matches(string(referenceModel.rxns), rxnsNames));

        samples = project.comparisons.(comparisonName).samplingComparison.orderedSamples;
        zeroRxns = find(sum(samples == 0, 2) == size(samples, 2));
        rxnIds = setdiff(rxnIds, zeroRxns);
        rxnsNames = referenceModel.rxns(rxnIds);
        
        orderedLb = getOrderedFeatureMatrix(project, reference, reference, "rxns", "model.lb");
        orderedUb = getOrderedFeatureMatrix(project, reference, reference, "rxns", "model.ub");
        orderedMappingRxnMatrix = getOrderedFeatureMatrix(project, modelNames, reference, "rxns", "mappedDiscretizedRxns");
    
        % check with samplings are not zero 
        if isfield(ref,'mediumComposition') %in case there is no medium defined
            rxnsNames = string(rxnsNames);
            % checking which column of the medium composition table entails
            % the rxn names -> just checking for the column with most
            % overlap
            isText = varfun(@(x) iscell(x) || isstring(x) || ischar(x), ...
                            ref.mediumComposition, ...
                            "OutputFormat", "uniform");

            textTable = ref.mediumComposition(:, isText);
            nMatches = varfun(@(x) sum(ismember(string(x), rxnsNames)), textTable);
            [maxMatches, bestColumn] = max(nMatches{:,:});
            bestColumnName = ref.mediumComposition.Properties.VariableNames{bestColumn};

            mediumConstrained = ismember(rxnsNames, ref.mediumComposition.(bestColumn)); %| ...
                             % ismember(rxns_names, ref.manual_set_boundaries.unwanted_export) | ...
                             % ismember(rxns_names, ref.manual_set_boundaries.unwanted_import);
        end

        orderedUb = orderedUb(rxnIds, :);
        orderedLb = orderedLb(rxnIds, :);
        orderedMappingRxnMatrix = orderedMappingRxnMatrix(rxnIds, :);
       
        rxnAbbr = referenceModel.rxns(rxnIds);
        rxnFormulas = string(printRxnFormula(referenceModel, rxnAbbr, false));
  
        % get rxn gene rules to add to the tabler
        symbolGPRrules = string(cellfun(@(rxnName)getRxnSymbolRule(project.models.(reference), rxnName), string(rxnsNames), 'UniformOutput', false));

        T = table(rxnFormulas, mediumConstrained,...
            join(string(orderedMappingRxnMatrix), "|", 2), symbolGPRrules, ...
                  'VariableNames', ["Reaction Formula", "Medium Constrained", ...
                  join(string(project.comparisons.(comparisonName).modelNames), "_"), ...
                  "Symbol GPR Rules"], 'RowNames', rxnsNames);
        
        T = T(flip(string(T.Properties.RowNames)), :);
        T.lb = flip(orderedLb);
        T.ub = flip(orderedUb);

        % Convert row names into a column
        T_display = T;
        T_display = addvars(T_display, string(T_display.Properties.RowNames), ...
                            'Before', 1, ...
                            'NewVariableNames', "Reaction");
        
        % Create UI figure
        fig = uifigure('Name', "Metabolite Fluxsum in subSystem: " + titleFig, ...
                       'Position', [100 100 1200 600]);
        
        % Create table filling the entire figure
        tbl = uitable(fig, ...
            'Data', T_display{:, :}, ...
            'ColumnName', T_display.Properties.VariableNames, ...
            'Units','normalized', ...
            'Position', [0 0 1 1], ...
            'FontSize', 18, ...
            'ColumnWidth', 'auto');
        end
        plotName = replace(plotName, ["_", "-", "/", " "], "");
        figs.("violinFlux" + plotName) = figStruct;

    end

end

function fluxCell = getFlux(project, comparisonName, rxnIdx,slot)
    
    arguments
        project
        comparisonName
        rxnIdx = []
        slot (1,1) string {mustBeMember(slot, ["orderedSamples","orderedllSamples"])} = "orderedSamples" 
    end
    
    reference = project.comparisons.(comparisonName).referenceModel;
    
    if isempty(rxnIdx)
        rxnIdx = {find(ones(length(project.models.(reference).model.rxns), 1))};
    end
    
    [samples, ~] = findFieldRecursive(project.comparisons.(comparisonName),slot);
    
    fluxCell = cellfun(@(x) samples(x, :), ...
                           rxnIdx, 'UniformOutput', false);

end


function datKLD_reordered = reorderKLDComparisons(groups, datKLD, datKLDlabels)

    nGroups = numel(groups);
    nComparisons = nchoosek(nGroups, 2);

    % Expected comparison order from the violin plot
    plotLabels = strings(nComparisons, 1);

    k = 1;
    for i = 1:nGroups-1
        for j = i+1:nGroups
            plotLabels(k) = groups(i) + "_vs_" + groups(j);
            k = k + 1;
        end
    end

    % Determine which input column belongs to each plot comparison
    columnOrder = zeros(1, nComparisons);

    for j = 1:nComparisons

        plotParts = split(plotLabels(j), "_vs_");

        for i = 1:numel(datKLDlabels)

            parts = split(datKLDlabels(i), "_vs_");

            % Match either direction
            if (parts(1) == plotParts(1) && parts(2) == plotParts(2)) || ...
               (parts(1) == plotParts(2) && parts(2) == plotParts(1))

                columnOrder(j) = i;
                break
            end
        end
    end

    % Reorder columns
    datKLD_reordered = datKLD(:, columnOrder);

end

function addSignificanceBars(dat,comparisonLabels, groups,dataViolin)

    % Data range
    data = dat;

    ymax = max(dataViolin, [], "all");
    ymin = min(dataViolin, [], "all");
    yrange = ymax - ymin;

    % Loop over comparisons
    nComparisons = numel(comparisonLabels);

    for k = 1:nComparisons

        % e.g. "Control_vs_StageI"
        comparison = comparisonLabels(k);

        % Split comparison into group names
        parts = split(comparison, "_vs_");

        group1 = parts(1);
        group2 = parts(2);

        % Find positions of groups on the violin plot
        pos1 = find(groups == group1);
        pos2 = find(groups == group2);

        % Make sure both groups exist
        if isempty(pos1) || isempty(pos2)
            error("Could not find groups for comparison: %s", comparison);
        end

        % Statistical test
        p = data(k);

        % Height of significance bar
        h = 0.02 * yrange;

        % Draw bracket
        plot([pos1 pos1 pos2 pos2], ...
             [ymax ymax-h ymax-h ymax], ...
             'k-', 'LineWidth', 1.5)

        % p-value
        text((pos1+pos2)/2, ymax+h+0.05*yrange, ...
             sprintf('p = %.3f', p), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom');


    end

    % Add space above the significance annotations
    currentYLim = ylim;
    
    ylim([ ...
        currentYLim(1), ...
        max(currentYLim(2), ymax + 0.15*yrange) ...
    ]);

end