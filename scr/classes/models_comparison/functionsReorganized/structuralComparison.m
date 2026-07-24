function structureAnalysis = structuralComparison(project, modelList, referenceModel)
    % The structure comparison is a function that compares the models
    % listed on structural differences. Structural differences in the
    % context of Fastcore can be defined as the set of reactions that are
    % kept when running fastcore. This means we check for the existence of
    % rxns, metabolites and genes in the model, and their overlap between
    % models. Reaction existence will be analysed in different sets.
    % Reaction exisitence in different metabolic subsystems/pathways,
    % existence of reaction depending on their initial classification of
    % expressed/not expressed core/not expressed. This figures are meant to
    % give an indication of which reactions were kept in the model, how
    % high the overlap of those are between the models, which of the
    % overlapping reactions are within the core or not expressed, which
    % pathways the genes in the outer and intersection are part of etc.
    % Input:
    %   - project: fastcore workflow project, with a run of 
    %              singleModelAnalysis already performed, and the active
    %              analysis already set using
    %              chooseActiveAnalysisForComparison
    %   - modelList: List of models to compare
    %   - referenceModel: for the structural comparison of wether a
    %                      reaction is existent or not, a reference model needs to be defined
    %                      therefore you need to define a reference model
    %                      which is also in the models slot of the project
    %                      object
    % Output:
    %   - structure analysis: yet to be structured properly #TODO
    arguments
        project (1,1) struct
        modelList (1,:) string
        referenceModel (1,1) string
    end

    modelsList = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    models = structfun(@(x) x.model, modelsList, 'UniformOutput', false);
    structureAnalysis.modelNames = string(fieldnames(models));


    % get model sizes - # genes,reactions and metabolites
    data = struct2array(structfun(@(x) {numel(x.rxns); numel(x.mets); numel(x.genes)}, ...
                           models, 'UniformOutput', false))';
    array2table(data, ...
                    'VariableNames', {'count_reactions','count_metabolites','count_genes'}, ...
                    'RowNames', string(fieldnames(models))')

    [rxnPresence, rxnMapping] = getOrderedFeatureMatrix(project,modelList, "rxns", referenceModel);
    structureAnalysis.rxn_mapping_table = array2table(rxnMapping, "VariableNames", modelList, "RowNames", string(project.models.(referenceModel).model.rxns));

    
    % -- Visualization: Discretization status for expression of genes in model on sample level, on model level as well as the mapping/discretization on rxn level
    % -> gives you a feeling of how many reactions in the model are from the core, how many of the rxns that were notExpressed made it in regardless etc.
    
    % check if the models has all the information needed to analyse the
    % core reactions
    
    
    if all(structfun(@(x) isfield(x, "discretized_data"), modelsList))
        % get the reaction mapping (sample and model level) as well as the discretization values for each reaction/gene in the model 
        replacementValue = "mappedDiscRxnsSample"; % get the fba solution values
        orderedMappingRxnMatrixSampleWise = int8(getOrderedFeatureMatrix(project, modelList, "rxns", referenceModel, replacementValue));
        replacementValue = "mappedDiscRxns"; % get the fba solution values
        orderedMappingRxnMatrix = int8(getOrderedFeatureMatrix(project, modelList, "rxns", referenceModel, replacementValue));
        replacementValue = "discretized_data.values"; % get the fba solution values
        orderedMappingExprDiscMatrix = int8(getOrderedFeatureMatrix(project, modelList, "genes", referenceModel, replacementValue));
        
        if all(structfun(@(x) isfield(x.settings, "script_parameters"), modelsList))
            if all(structfun(@(x) isfield(x, "sample_metadata"), modelsList))
                if all(structfun(@(x) isfield(x.settings.script_parameters, "columns_to_define_model_samples_on"), modelsList))
                    % get the names of the single samples from the metadata slot - used in the following plots
                    columnNames = struct2cell(structfun(@(x)  string(x.sample_metadata{:,1}) + "_" + ...
                                                  string(x.sample_metadata.(x.settings.script_parameters.columns_to_define_model_samples_on)),...
                                            modelsList, "UniformOutput", false));
                    columnNames = vertcat(columnNames{:});
                else
                    columnNames = string(1:size(orderedMappingRxnMatrixSampleWise,2));
                end
            else
                columnNames = string(1:size(orderedMappingRxnMatrixSampleWise,2));
            end
        else
            columnNames = string(1:size(orderedMappingRxnMatrixSampleWise,2));
        end
        
        % get the data into one object to loop over
        datasets = {orderedMappingExprDiscMatrix, orderedMappingRxnMatrixSampleWise, orderedMappingRxnMatrix};   % replace with your actual dataset variables
        datasetNames = ["ordered_mapping_rxn_matrix_sample_wise", "ordered_mapping_expr_disc_matrix", "ordered_mapping_rxn_matrix"];  % optional titles
        xlabelsPlots = ["Samples", "Samples", "Models"]; 
        xticksPlots = {columnNames, columnNames, string(fieldnames(modelsList))}; 
        ylabelsPlots = ["# genes for genes which are in the context specific models", "# reactions for reactions which are in the context specific models", "# reactions for reactions which are in the context specific models"];  
        titlePlots = ["after discretization: ", "after mapping the gpr rules: ", "after applying consensus proportion"];  
    
        
    
        % Determine all unique discretization values across datasets (excluding 13)
        % the value 13 has been set to indicate that the rxn/gene is not in the
        % model, so the discretization is not shown, in these figures only the
        % discretization is shown of the genes/rxns in the model, the figures
        % for all genes, rxns are done in the QC script when preparing the data
        % for the model creation !!
        plots.dataDiscretization = getDiscretizationHist(datasets);
    end     
    if all(structfun(@(x) isfield(x, "core_reactions"), modelsList))

        %%% Visualize the core reaction per model
        data = struct2cell(structfun(@(x) [length(x.core_reactions) - sum(ismember(x.core_reactions, x.model.rxns)); ...
                                                sum(ismember(x.core_reactions, x.model.rxns));...
                                                length(x.model.rxns) - sum(ismember(x.core_reactions, x.model.rxns));...
                                                length(x.model.rxns)], ...
                                                modelsList, 'UniformOutput', false));
        data = [data{:}];
    
        plots.coreReactions = plotCoreInModel(data, modelsList);
    
        % -- Visualization: Looking in deeper into the core reactions, the core
        % is what is defined by the data, therefore portrays the underlying
        % biological chnages, so the question is which reactions are part of
        % the outer and intersections we saw in the previous venn/intersection
        % diagramm ? are the differences in core reactions only due to
        % exchange/import ? transporters ? This should be avoided!
    
        % create an upsetr plot for the all the inter and outersections
        % filter out the main intersection -> the one with the longest name
    
        plots.coreReactionsIntersections = getUpsetPlotCore(project, modelsList, structureAnalysis);
        % this function only works with a comparison of up to 4 models!!
        % otherwise the plots get too complex!

    end


    % -- Visualization: Get the jaccard similarity on basis of the
    % gene,metabolite and reaction presence in the corresponding models
    % How similar are my models structuraly, which models are more similar
    % to each other than others ?

    for fieldToInvestigate = ["genes", "mets", "rxns"]
        [orderedFeature, ~] = getOrderedFeatureMatrix( ...
            project, modelList, fieldToInvestigate, referenceModel);

        % Plot Venn / Heatmap of intersections based on presence
        plots.intersections.(fieldToInvestigate) =  plotFlexibleVenn( ...
                                                                    orderedFeature, ...
                                                                    structureAnalysis.modelNames, ...
                                                                    "Structural model comparison: " + fieldToInvestigate + " presence",...
                                                                    "visiblePlot", "off");

        % get the jaccard distances - based on reaction presence
        % Compute Jaccard distances

        plots.jaccardDist.(fieldToInvestigate) =  plotJaccard( ...
                                                                 orderedFeature, ...
                                                                 structureAnalysis.modelNames, ...
                                                                 "Jaccard similarity of " + fieldToInvestigate + " presence (0 or 1) between models",...
                                                                 "visiblePlot", "off");





    end

    % -- Visualization: Get reaction presence for each model in comparison
    % to the defined reference model -> visualization per subsystem
    % Where does the difference I see in the jaccard plot come from ? form
    % which subsystem, which subsystem is most different in pairwise
    % comparison ? 

    plots.reactionPathwayPresence = pathwayPresenceHeat(project, referenceModel);

    structureAnalysis.plots.struct = plots;


    function fig = getDiscretizationHist(datasets)
        % This function produces a nested barplot that shows how many of
        % rxns have a discretization status of 0 1 or -1.For the different
        % models but also for the different samples in the models!
        numDatasets = length(datasets);
        allValues = [];
        for k = 1:numDatasets
            allValues = union(allValues, setdiff(unique(datasets{k}), 13));
        end
        allValues = sort(allValues);  % e.g., [-1 0 1]
        
        % Define a colormap with one color per value
        cmap = lines(length(allValues));
        
        % Create figure
        fig = figure('Color','w','Position', [100 100 2000*numDatasets 2000],...
                                           'Visible', 'off');
        
        for k = 1:numDatasets
            dataset = datasets{k};
            xlabelPlot = xlabelsPlots(k);
            xtickPlot = xticksPlots{k};
            ylabelPlot = ylabelsPlots(k);
            titlePlot = titlePlots(k);
            
            % counts per category per sample (make sure all_values are included)
            counts = zeros(length(allValues), size(dataset,2));
            for i = 1:length(allValues)
                counts(i,:) = sum(dataset == allValues(i), 1);
            end
        
            % stacked barplot in subplot
            ax = subplot(1, numDatasets, k);
            b = bar(ax, counts', 'stacked');
        
            % Assign consistent colors
            for i = 1:length(allValues)
                b(i).FaceColor = cmap(i, :);
            end
        
            % percentages
            tot = sum(counts, 1);
            pct = 100 * counts ./ tot;
            pct = round(pct);
        
            % write percentages inside bars
            for i = 1:size(counts, 2)
                y0 = 0;
                for j = 1:size(counts, 1)
                    if counts(j, i) > 0
                        text(i, y0 + counts(j, i)/2, ...
                            sprintf('%g%%', pct(j, i)), ...
                            'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'middle', ...
                            'FontSize', 13, 'Color', 'w', 'FontWeight', 'bold');
                    end
                    y0 = y0 + counts(j, i);
                end
            end
        
            % axes labels and formatting
            ax.FontSize = 14;
            xlabel(xlabelPlot, 'FontSize', 16)
            ylabel(ylabelPlot, 'FontSize', 16)
            title(titlePlot, 'FontSize', 18)
        
            xticks(1:length(xtickPlot))
            xticklabels(regexprep(xtickPlot, "_", "-"))
        end
        
        % Single legend for the whole figure
        lgd = legend(string(allValues), 'Location', 'northeastoutside');
        lgd.FontSize = 14;
        lgd.Title.String = "Discretization status";
    end
    
    function plt = pathwayPresenceHeat(project, referenceModel)
        % This function produces a heatmap showing the rxns presence in
        % comparison to a reference model in the color scale (0-1) and the
        % absolute number of rxns per pathway in every model as text in the
        % tiles, additionally it displays the pathway size of each pathway
        % on the left side in an additional barplot!
        pathways = string(project.models.(referenceModel).model.subSystems); % get pathways from reference model
        uniquePathways = unique(pathways); 
    
        % for every pathway get the matrix identifying the rnxs from reference
        % model in this pathway
        groups = arrayfun(@(x) find(pathways == x), uniquePathways, 'UniformOutput', false);
        numGroups = length(groups);
        G = zeros(numGroups, size(orderedFeature,1));
    
        for g = 1:numGroups
            G(g, groups{g}) = 1;
        end
    
        % get the presence per subsystem in the context specific models 
        % by mulitplying the rxns presence for each subsystem (matrix ordered features) 
        % with the matrices defining the subsystem for every rxns
        
        pathwayCounts = array2table(G * orderedFeature, ...
                                     'VariableNames', structureAnalysis.modelNames,...
                                     'RowNames', string(cellstr(uniquePathways)));
        % add reference model count to be able to make a relative abundance
        pathwayCounts.referenceModel = groupcounts(pathways);
    
        relativeCounts = array2table(pathwayCounts{:, 1:end-1} ./ pathwayCounts.referenceModel, ...
                                     'VariableNames', structureAnalysis.modelNames,...
                                     'RowNames', cellstr(uniquePathways));
    
        % get the idx of the most variant pathways in terms of rxns presence
        relativeCounts.row_var = var(relativeCounts{:,:}, 0, 2);
        pathwayCounts.row_var = var(pathwayCounts{:, 1:end-1}, 0, 2);
        pathwayCounts = pathwayCounts(pathwayCounts.referenceModel < 1000,:);
        % Get indices of top n highest variance rows
    
    
        pathwayCounts = sortrows(pathwayCounts,"row_var","descend");
        pathwayCounts = pathwayCounts(find(pathwayCounts.row_var ~= 0),:);
        relativeCounts = relativeCounts(pathwayCounts.Properties.RowNames,:);
        % plot top 20 most variant pathways between the choosen models
        data = relativeCounts{:, 1:end-1};
        rowNames = string(relativeCounts.Properties.RowNames);
        colNames = structureAnalysis.modelNames;
    
        %%%%%%%%%%
    
        plt = figure('Color', 'w', 'Position', [20 20 700 300], 'Visible', 'off');
        tiledlayout(1, 4, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact')
    
        % ---- Bar plot (LEFT) ----
        ax1 = nexttile(1);
        barh(pathwayCounts.referenceModel)
        title('Subsystem size in the reference model')
        xlabel('# rxns in the reference model')
    
        ax1.YTick = 1:numel(rowNames);
        ax1.YTickLabel = rowNames;
        ax1.YDir = 'reverse';
        ax1.YAxisLocation = 'right';   % ⭐ labels between plots
        ax1.TickLength = [0 0];        % removes tick marks
        %ax1.XTick = [];               % remove x ticks
        ax1.YTick = 1:numel(rowNames);
        ax1.YTickLabel = rowNames;    % keep labels
        ax1.YAxisLocation = 'right';  % labels between plots
        ax1.FontSize = 12;   % bar plot labels
        % Flip the horizontal axis
        ax1.XDir = 'reverse';
    
    
    
    
        % ---- Heatmap (RIGHT, spanning 2 tiles) ----
        % z-scaling of the data -> so that the colorod
        ax2 = nexttile(2, [1 3]);
    
        imagesc(data)
        cmap = getColorPallette();
        colorbar
        title("Relative counts of subsystem rxn occurence/reference model")    % grayscale
        ax2.XTick = 1:numel(colNames);
        ax2.XTickLabel = colNames;
    
        ax2.YTick = 1:numel(rowNames);
        ax2.YTickLabel = [];    % labels shown only once
    
        ax2.TickLength = [0 0];
    
        xlabel('Models')
        ax2.FontSize = 12;   % heatmap labels
    
        % Get size of the data
        [nRows, nCols] = size(data);
    
    
        % Loop over every cell and place the absolute number from pathway_counts
        for i = 1:nRows
            for j = 1:nCols
                % You want absolute numbers, not relative counts, so use pathway_counts
                % (or multiply relative_counts by referenceModel if needed)
                value = pathwayCounts{i, j}; % +1 because first column is referenceModel
                % Place text at the center of the tile
                text(ax2, j, i, num2str(value), ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'Color','k', ...          % black text
                    'FontSize',10)
            end
        end
        % ---- Align rows ----
        linkaxes([ax1 ax2], 'y')

    end
    
    function fig = plotCoreInModel(data, modelsList)
        % This function visualizes how many of the core rxns defined by
        % rFastcormics make it in and how many of the core rxns are
        % overlapping between the models.=
        upperData = data([3, 2], :);
        categories = fieldnames(modelsList)'; 
    
        % -- Visualization: get a sense of how many of the core reactions made
        % it into your model. In theory 100 percent of the reactions should be
        % in, but in practice this is not gona happen, adjusting the
        % thresholding/discretization in the preprocessing of the data for
        % rfastcormics could change the precentages you see in the following
        % figure
    
    
    
        fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1500 1500]);
        tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
    
        % --- first barplot
        ax1 = nexttile(1);
        bar(upperData', 'stacked')
    
    
        % Labels
        set(gca, 'XTickLabel', categories, 'FontSize', 14)  % increase tick labels font
        xlabel('Models', 'FontSize', 14)
        ylabel('# rxns', 'FontSize', 14)
    
        % Legend
        legend({"Non-core reactions","Core reactions"}, 'Location', 'northwest', 'FontSize', 14)
    
        % Title
        title('Core and non-core reactions per model', 'FontSize', 14)
    
        % ---- Compute percentages of upper stack ----
        total = sum(upperData, 1);                   % total per model
        percentUpper = 100 * upperData(2, :) ./ total;  % percentage of upper bar
    
        % ---- Add text labels on top of upper bars ----
        for i = 1:size(upperData, 2)  % loop over models
            % x-position is the bar center, y-position is height of lower + upper
            x = i;
            y = upperData(1, i) + upperData(2, i)/2;  % middle of upper stack
            text(x, y, sprintf('%.1f%%', percentUpper(i)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold');
        end
    
         % --- second barplot
         data = data(1:2, :);  % only core vs non-core counts
        ax2 = nexttile(2);
        hb = bar(data', 'stacked');
    
        % Labels
        set(gca, 'XTickLabel', categories, 'FontSize', 14)
        xlabel('Models', 'FontSize', 14)
        ylabel('# rxns', 'FontSize', 14)
    
        % Legend
        legend({ "Not included","Included"}, 'Location', 'northwest', 'FontSize', 14)
    
        % Title
        title('Core reactions used into model building.', 'FontSize', 14)
    
        % ---- Compute percentages of upper stack ----
        total = sum(data, 1);                   % total per model
        percentUpper = 100 * data(2,:) ./ total;  % percentage of upper bar
    
        % ---- Add text labels on top of upper bars ----
        for i = 1:size(data, 2)  % loop over models
            % x-position is the bar center, y-position is height of lower + upper
            x = i;
            y = data(1, i) + data(2, i)/2;  % middle of upper stack
            text(x, y, sprintf('%.1f%%', percentUpper(i)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold');
        end
    
        % Tile 3 (THIS IS THE KEY)
        ax3 = nexttile(3, [1 2]);   % column 2, span both rows
        axis(ax3, 'off')
        hold(ax3, 'on')
    
    
        core_reactionsIncluded = struct2cell(structfun(@(x) x.core_reactions(find(ismember(x.core_reactions, x.model.rxns)))', ...
                                                modelsList, 'UniformOutput', false));
        core_reactionsIncluded = unique([core_reactionsIncluded{:}]);
    
        corePresence = structureAnalysis.rxn_mapping_table{core_reactionsIncluded,:} ~= 0;
        [figV,idxInterOutersections,~] = plotFlexibleVenn(corePresence, structureAnalysis.modelNames, ... 
                                                            "Structural model comparison: core rxns presence", "visiblePlot", "off");
    
    
        if string(class(figV)) == 'matlab.ui.Figure'
    
            % Find axes inside Venn figure
            axV = findobj(figV,'Type','axes','-not','Tag','legend');
            axV = axV(1);
    
            % Copy graphics
            copyobj(allchild(axV), ax3)
    
            % Fix geometry
            axis(ax3, 'tight')
            axis(ax3, 'equal')
            ax3.Clipping = 'off';
            close(figV)
        else 
            % store main figure handle
            mainFig = gcf;
    
            % create Venn/heatmap figure
            [figV, idxInterOutersections,~] = plotFlexibleVenn( ...
                corePresence, structureAnalysis.modelNames, ...
                "Structural model comparison: core rxns presence");
    
            % extract data
            X = figV.XData;
            Y = figV.YData;
            C = figV.ColorData;
    
            % close the temporary figure
            close(ancestor(figV, 'figure'))
    
            % activate main figure again
            figure(mainFig)
    
            % place heatmap in tile
            nexttile(3, [1 2])
            heatmap(X, Y, C)
    
            title('Structural model comparison: core reactions presence')
        end
    end

    function plt = getUpsetPlotCore(project, modelsList, structureAnalysis)
        % This figure shows the pathway of the inter and outersections
        % between the core reactions in the model. This is a quality figure
        % to controll in which subsystems we see a difference based on the
        % defined core of the different models. We want more than just the
        % Transport or exchangers to be different between models, in this
        % figure we can control that!
        % The figure only works for comparisons of 4 or less models, since
        % the Venn diagramm calculating the intersections and outersections
        % only works for 4 or less, otherwise the histogramm would get to
        % big to visualize. 
        
        core_reactionsIncluded = struct2cell(structfun(@(x) x.core_reactions(find(ismember(x.core_reactions, x.model.rxns)))', ...
                                                modelsList, 'UniformOutput', false));
        core_reactionsIncluded = unique([core_reactionsIncluded{:}]);
        corePresence = structureAnalysis.rxn_mapping_table{core_reactionsIncluded, :} ~= 0;
        [figV,idxInterOutersections,~] = plotFlexibleVenn(corePresence, structureAnalysis.modelNames, ... 
                                                            "Structural model comparison: core reactions presence", "visiblePlot", "off");

        if string(class(figV)) == 'matlab.ui.Figure'

            namesIntersections = fieldnames(idxInterOutersections);
            [~, allIntersection] = max(cellfun(@(x) length(x), namesIntersections));
            idxInterOutersections = rmfield(idxInterOutersections, namesIntersections(allIntersection));
            % now get the pathway of every entry
            interOutersectionsPathways = structfun(@(x) string(project.models.(referenceModel).model.subSystems(x)),...
                                                     idxInterOutersections, 'UniformOutput', false);
            C = struct2cell(interOutersectionsPathways);
            uniquePathways = unique(vertcat(C{:}));
    
    
            % Preprocess pathways: collapse transport
            S = structfun(@(x) regexprep(x, "^Transport.*", "Transport"), interOutersectionsPathways, 'UniformOutput', false);
            pathwaysUnique = unique(regexprep(uniquePathways, "^Transport.*", "Transport"));
    
            barNames = string(fieldnames(S));
            nBars = numel(barNames);
    
            % Build count matrix
            Y = cellfun(@(b) sum(S.(b)' == pathwaysUnique, 2)', barNames', 'UniformOutput', false);
            Y = cat(1, Y{:});
    
            % Sort bars by total counts (descending)
            [~, sortIdx] = sort(sum(Y, 2), 'descend');
            Y = Y(sortIdx, :);
            barNamesSorted = barNames(sortIdx);
    
            % Plot
            plt = figure('Color', 'w', 'Position', [100 100 6000 2000], 'Visible', 'off');
            b = bar(Y, 'stacked');
    
            % Generate a qualitative colormap with enough colors
            numColors = size(Y, 2);
            % Example 20-color qualitative palette (from ColorBrewer / Tableau)
            cmap = [ ...
                166 206 227;
                31 120 180;
                178 223 138;
                51 160 44;
                251 154 153;
                227 26 28;
                253 191 111;
                255 127 0;
                202 178 214;
                106 61 154;
                255 255 153;
                177 89 40;
                141 211 199;
                255 255 179;
                190 186 218;
                251 128 114;
                128 177 211;
                253 180 98;
                179 222 105;
                252 205 229] / 255;  % Normalize 0-1
    
            % Apply colors to each category
            for k = 1:numColors
                b(k).FaceColor = cmap(mod(k-1, size(cmap, 1)) + 1, :);
            end
    
    
            % Labels and legend
            ax = gca;
            ax.XTickLabel = regexprep(barNamesSorted, "_", " ");
            ax.FontSize = 20;
            xlabel('Model intersections/outersections', 'FontSize', 20)
            ylabel('# Core Reactions', 'FontSize', 20)
            title('Count of Core reactions per pathway and intersection/outersection', 'FontSize', 20)
    
            lgd = legend(pathwaysUnique, 'Location', 'northeast');
            lgd.FontSize = 20;
        else
            plt = [];
        end
    end
end