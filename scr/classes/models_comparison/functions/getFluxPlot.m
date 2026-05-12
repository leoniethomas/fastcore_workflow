function fig = getFluxPlot(project,comparison_name, idxToVis,options)
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
    %                           - reducedCost = true (default false) to display the reduced
    %                             cost values for every reaction as the color
    %                             of the FBA dot.
    %                           - thresholdFlux= wether to apply an upper
    %                             lower ,no upper or lower or to include 
    %                             also the fba =0 reactions(all) to the selected
    %                             reaction fba values
    %                           - 
    % 
    % Output:                   Display of figure                
    % 

    arguments
        project
        comparison_name (1,1) string
        idxToVis
        options.FVA  (1,1) logical = false
        options.reducedCost (1,1) logical = false
        options.thresholdFlux (1,1) string {mustBeMember(options.thresholdFlux,["lower", "upper","none","all"])} ="none" 
        options.titlePlots = ""
        options.visiblePlots ="on"
    end


    modelList = project.comparisons.(comparison_name).modelNames;
    referenceModel = project.comparisons.(comparison_name).referenceModel;
    
    replacement_value = "analysis.active.FBA.v"; % get the fba solution values
    ordered_fba_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacement_value);
    replacement_value = "mappedDiscRxns"; % get the fba solution values
    ordered_mapping_rxn_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacement_value);
    
    
    if options.thresholdFlux == "upper"
        get_exchange_rxns_idx = intersect(find(sum(ordered_fba_matrix,2) ~=0 & sum(ordered_fba_matrix < 0,2) ~= 0), ...
                                          idxToVis);  
        title_word = "negative flux reactions";
    elseif options.thresholdFlux == "lower"
        get_exchange_rxns_idx = intersect(find(sum(ordered_fba_matrix,2) ~=0 & sum(ordered_fba_matrix > 0,2) ~= 0), ...
                                                      idxToVis); 
        title_word = "positive flux reactions";
    elseif options.thresholdFlux == "none"
        get_exchange_rxns_idx = intersect(find(sum(ordered_fba_matrix,2) ~=0), ...
                                                      idxToVis);

        title_word = options.titlePlots;
    elseif options.thresholdFlux == "all"
        get_exchange_rxns_idx = idxToVis;
        title_word = options.titlePlots;
    else
        error("wrong value for threshold choosen. Possible values: lower, upper, none")
    end

    ordered_fba_matrix_ex = ordered_fba_matrix(get_exchange_rxns_idx,:);
    ordered_mapping_rxn_matrix_ex = ordered_mapping_rxn_matrix(get_exchange_rxns_idx,:);
    rxn_names = project.models.(referenceModel).model.rxns(get_exchange_rxns_idx);

    if options.FVA
        replacement_value = "analysis.active.FVA.maxFlux"; % get the fba solution values
        ordered_fva_max_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacement_value);
        replacement_value = "analysis.active.FVA.minFlux"; % get the fba solution values
        ordered_fva_min_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacement_value);
        ordered_fva_max_matrix_ex = ordered_fva_max_matrix(get_exchange_rxns_idx,:);
        ordered_fva_min_matrix_ex = ordered_fva_min_matrix(get_exchange_rxns_idx,:);
        
        if options.thresholdFlux == "all"
            idx_FVA_var = ~(all(ordered_fva_max_matrix_ex == 0,2) & all(ordered_fva_max_matrix_ex == 0,2));
            ordered_fva_max_matrix_ex = ordered_fva_max_matrix_ex(idx_FVA_var,:);
            ordered_fva_min_matrix_ex = ordered_fva_min_matrix_ex(idx_FVA_var,:);
            rxn_names = rxn_names(idx_FVA_var);
            ordered_fba_matrix_ex = ordered_fba_matrix_ex(idx_FVA_var,:);
            ordered_mapping_rxn_matrix_ex = ordered_mapping_rxn_matrix_ex(idx_FVA_var,:);
        
        end
        
        if options.reducedCost

            models = rmfield(project.models,setdiff(fieldnames(project.models), modelList));
            if all(structfun(@(x) isfield(x.analysis.active.FBA,"w"),models)) && all(structfun(@(x) ismember('one',x.analysis.active.parameters.Value(x.analysis.active.parameters.Analysis == "FBA")), models))
                
                for m = modelList'
                    nRxns = length(project.models.(m).model.rxns);
                    w     = project.models.(m).analysis.active.FBA.w;  % 7350×1
                    project.models.(m).analysis.active.FBA.reducedcost = w(2*nRxns+1 : 3*nRxns);    % net flux variables ← use this one
                end
                replacement_value = "analysis.active.FBA.reducedcost"; % get the fba solution values
                ordered_reducedCost_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacement_value);
                ordered_reducedCost_matrix_ex = ordered_reducedCost_matrix(get_exchange_rxns_idx,:);
                if options.thresholdFlux == "all"
                    ordered_reducedCost_matrix_ex = ordered_reducedCost_matrix_ex(idx_FVA_var,:);
                end
            else
                error("The reduced costs could not be found in the .w slot of the FBA analysis. In case you did not use the 'one' minNorm for the FBA the shadowprices might be stored elsewere!")
            end
        end

        
    end




    
    rxn_formulas = printRxnFormula(project.models.(referenceModel).model, 'rxnAbbrList', rxn_names, 'printFlag', false);

    % medium composition 
    % do the models have the same  medium ? 
    % if so I can add it as column to the plot, otherwise it is not added
    models = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    media_for_models = structfun(@(x) x.settings.medium, models);
    ref = fieldnames(models);
    ref = models.(ref{1,1}).settings.medium;
    medium_is_equal_between_models = all(arrayfun(@(x) isequaln(x, ref), media_for_models));
    

    if medium_is_equal_between_models
        medium_constrained = ismember(rxn_names, ref.medium_composition.ExRxns_Recon3D);% | ...
                         %ismember(rxn_names, ref.manual_set_boundaries.unwanted_export) | ...
                         %ismember(rxn_names, ref.manual_set_boundaries.unwanted_import);

        

        ordered_lb = getOrderedFeatureMatrix(project,referenceModel,...
                                             "rxns",referenceModel,"model.lb");
        ordered_ub = getOrderedFeatureMatrix(project,referenceModel,...
                                             "rxns",referenceModel,"model.ub");
        ordered_ub = ordered_ub(get_exchange_rxns_idx,:);
        ordered_lb = ordered_lb(get_exchange_rxns_idx,:);
        if options.FVA & options.thresholdFlux == "all"
            ordered_lb = ordered_lb(idx_FVA_var,:);
            ordered_ub = ordered_ub(idx_FVA_var,:);
        end

        % get rxn gene rules to add to the table
        symbol_gpr_rules = string(cellfun(@(rxnName)getRxnSymbolRule(project.models.(referenceModel),...
                                                       rxnName),string(rxn_names),'UniformOutput', false));
       
        T = table(rxn_formulas, medium_constrained,...
            join(string(ordered_mapping_rxn_matrix_ex), "|", 2),symbol_gpr_rules,...
                  'VariableNames', ["Reaction Formula","medium constrained",...
                                    join(string(project.comparisons.(comparison_name).modelNames),"_"),...
                                    "symbol gpr rules"], ...
                  'RowNames',rxn_names);
        T = T(flip(string(T.Properties.RowNames)),:);
        T.lb = flip(ordered_lb);
        T.ub = flip(ordered_ub);
    end
    


    
    
    %if options.shadowPrice
    %    replacement_value = "analysis.active.FBA.basis.dual"; % get the fba solution values
    %    % shadow prices are measured for every metabolite therefore mapped according to the mets field
    %    ordered_shadowPrices_matrix = getOrderedFeatureMatrix(project,modelList,"mets",referenceModel,replacement_value);
    %end

    %%% Parameters for the figure 

    %%% =========================
    % Threshold & masks
    % =========================
    mu = mean(ordered_fba_matrix_ex, 2);
    
    x  = abs(mu);
    
    % 1D k-means to minimize within-group distances

    assert(length(x) >1,...
           'In the end only one rxn value was returned after filtering for the rxns that are active in your FBA solution, in order to visualize this plot, you need to formulate more rxns in rxnID or apply other threholds (try "all" threshold)!')
    
    [idx, C] = kmeans(x, 2, 'Replicates', 10);
    
    % Identify low-flux (dense) cluster
    [~, lowCluster] = min(C);
    
    isLow  = idx == lowCluster;
    isHigh = ~isLow;
    
    % Derive lowRange for documentation
    %lowRange = [-2 2];
    lowRange = [-max(x(isLow)), max(x(isLow))];

    
    isLow  = mu > lowRange(1) & mu < lowRange(2);
    isHigh = ~isLow;
    
    dataLow  = ordered_fba_matrix_ex(isLow , :);
    dataHigh = ordered_fba_matrix_ex(isHigh, :);
    
    rxn_names_low  = rxn_names(isLow);
    rxn_names_high = rxn_names(isHigh);
    
    %%% =========================
    % Adaptive relative heights
    % =========================
    alpha = 0.75;
    nLow  = numel(rxn_names_low);
    nHigh = numel(rxn_names_high);
    
    w = [nLow nHigh].^alpha;
    usableHeight = 1 - 0.08 - 0.10 - 0.03;
    
    heightLow  = usableHeight * w(1) / sum(w);
    heightHigh = usableHeight * w(2) / sum(w);
    
    bottomLow  = 0.10;
    bottomHigh = bottomLow + heightLow + 0.03;
    %%%
    

    fig = uifigure('Name', title_word + " with Table", ...
                       'Position',[100 100 1000 450],'Visible',options.visiblePlots);
        
    plotWidth = 0.52;
    
    %%% =========================
    % Helper for axes formatting
    % =========================
    fmtAxes = @(ax) set(ax, ...
        'FontSize',16, ...
        'PositionConstraint','outerposition', ...
        'GridColor',[.8 .8 .8], ...
        'GridAlpha',.5);
    
    if ~options.FVA & ~options.reducedCost % specify which of the fields need to be true and false!!!
        % in the case that only the FBA solution should be visualized we
        % use a grouped horizontal barplot to do so

        
        
        %%% =========================
        % TOP AXIS — high values
        % =========================
        if heightHigh ~= 0
            axTop = uiaxes(fig,'Units','normalized', ...
                'Position',[0.02 bottomHigh plotWidth heightHigh]);
            fmtAxes(axTop)

            if size(dataHigh, 1) == 1
                dataHigh = [dataHigh; NaN(1, size(dataHigh, 2))];
                dataHighOneBar = 1;
            else
                dataHighOneBar = 0;
            end
            
            barh(axTop, dataHigh, 'grouped');
            grid(axTop,'on')
            
            %axTop.YTick = [];
            %axTop.YColor = 'none';
            
            yticks(axTop, 1:numel(rxn_names_high))
            yticklabels(axTop, strrep(rxn_names_high, "_", "\_"))
            title(axTop, title_word)
            if heightLow == 0
                xlabel(axTop, 'Flux value [mMol/(gDW*h)]')
            end
            if dataHighOneBar
                ylim(axTop, [0.5, size(dataHigh, 1) - 0.4 ])
            end
        end

        %%% =========================
        % BOTTOM AXIS — low values
        % =========================
        if heightLow ~= 0
            axBottom = uiaxes(fig,'Units','normalized', ...
                'Position',[0.02 bottomLow plotWidth heightLow]);
            fmtAxes(axBottom)
            
            if size(dataLow, 1) == 1
                dataLow = [dataLow; NaN(1, size(dataLow, 2))];
                dataLowOneBar = 1; 
            else
                dataLowOneBar = 0;
            end
            barh(axBottom, dataLow, 'grouped');
            grid(axBottom,'on')
            
            yticks(axBottom, 1:numel(rxn_names_low))
            yticklabels(axBottom, strrep(rxn_names_low, "_", "\_"))
            if heightHigh == 0
                title(axBottom, title_word)
            end
            xlabel(axBottom, 'Flux value [mMol/(gDW*h)]')
            if dataLowOneBar
                ylim(axBottom, [0.5, size(dataLow, 1) - 0.4 ])
            end
        end
        
        %%% =========================
        % Reverse direction if needed
        % =========================
        if options.thresholdFlux == "upper"
            set([axTop axBottom], 'XDir','reverse')
        end
        
        %%% =========================
        % Legend
        if nHigh < nLow
            legend(axBottom, modelList, 'Location','best');
        else
            legend(axTop, modelList, 'Location','best');
        end
        
        
        


    else
        Q1_low  = ordered_fva_min_matrix_ex(isLow , :);
        MED_low = ordered_fba_matrix_ex(isLow , :);
        Q3_low  = ordered_fva_max_matrix_ex(isLow , :);
        
        Q1_high  = ordered_fva_min_matrix_ex(isHigh , :);
        MED_high = ordered_fba_matrix_ex(isHigh , :);
        Q3_high  = ordered_fva_max_matrix_ex(isHigh , :);
            
        % =========================
        % Check if reduced cost coloring is enabled
        % =========================
        add_color_legend = options.reducedCost;
        
        if add_color_legend
            cmap = colormap('cool');  % N colors
            N = size(cmap,1);
        
            % Split reduced cost into high and low
            reducedCost_low  = ordered_reducedCost_matrix_ex(isLow, :);
            reducedCost_high = ordered_reducedCost_matrix_ex(isHigh, :);
        
            % Helper function to scale matrix to colormap indices
            scaleToCmap = @(mat, valMin, valMax) round( (mat - valMin) / (valMax - valMin) * (N-1) ) + 1;
            
            % Get global min/max from full matrix
            valMin = min(ordered_reducedCost_matrix_ex(:));
            valMax = max(ordered_reducedCost_matrix_ex(:));
            
            if valMin == valMax
                % Entire matrix has one unique value — all same color
                scaledIdx_low  = ones(size(reducedCost_low));
                scaledIdx_high = ones(size(reducedCost_high));
            else
                % Check each subset individually
                if length(unique(reducedCost_low(:))) == 1
                    scaledIdx_low = ones(size(reducedCost_low));
                else
                    scaledIdx_low = scaleToCmap(reducedCost_low, valMin, valMax);
                end
            
                if length(unique(reducedCost_high(:))) == 1
                    scaledIdx_high = ones(size(reducedCost_high));
                else
                    scaledIdx_high = scaleToCmap(reducedCost_high, valMin, valMax);
                end
            end
        else
            add_color_legend = 0;
        end

       

       %%
        
        % --- TOP AXES: high flux ---
        %%% =========================
        % Plot high reactions
        %%% =========================
        if nHigh > 0
            axTop = uiaxes(fig,'Units','normalized','Position',[0.02 bottomHigh plotWidth heightHigh]);
            hold(axTop,'on'); fmtAxes(axTop);

            % --- Grid & font ---
            grid(axTop,'on')
            axTop.GridColor = [0.8 0.8 0.8];
            axTop.GridAlpha = 0.5;
            axTop.FontSize = 16;
        
            nGroups = nHigh;
            nPerGroup = size(Q1_high,2);
            boxHeight = 0.2; groupSep = 1;
            greyLevels = linspace(0.9,0.6,nPerGroup);
            greyRGBs = [greyLevels', greyLevels', greyLevels'];
        
            for i = 1:nGroups
                yBase = i * groupSep;
                for j = 1:nPerGroup
                    offset = (j-(nPerGroup+1)/2)*(boxHeight*1.3);
                    rectangle(axTop, 'Position',[Q1_high(i,j), yBase+offset-boxHeight/2, Q3_high(i,j)-Q1_high(i,j), boxHeight],...
                              'FaceColor',greyRGBs(j,:),'EdgeColor','none');
                    % Median dot
                    if add_color_legend
                        plot(axTop, MED_high(i,j), yBase+offset, 'o', ...
                        'MarkerSize',5, ...
                        'MarkerFaceColor', cmap(scaledIdx_high(i,j),:), ...
                        'MarkerEdgeColor', cmap(scaledIdx_high(i,j),:));

                    else
                        plot(axTop, MED_high(i,j), yBase+offset,'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
                    end
                end
            end
        
            yticks(axTop, 1:nGroups*groupSep)
            yticklabels(axTop, strrep(rxn_names_high,"_","\_"))
            title(axTop, title_word)
        end
        
        %%% =========================
        % Plot low reactions
        %%% =========================
        if nLow > 0
            
            axBottom = uiaxes(fig,'Units','normalized','Position',[0.02 bottomLow plotWidth heightLow]);
            hold(axBottom,'on'); fmtAxes(axBottom);
            
            % --- Grid & font ---
            grid(axBottom,'on')
            axBottom.GridColor = [0.8 0.8 0.8];
            axBottom.GridAlpha = 0.5;
            axBottom.FontSize = 16;
        
            nGroups = nLow;
            nPerGroup = size(Q1_low,2);
            boxHeight = 0.2; groupSep = 1;
            greyLevels = linspace(0.9,0.6,nPerGroup);
            greyRGBs = [greyLevels', greyLevels', greyLevels'];
            
            for i = 1:nGroups
                yBase = i * groupSep;
                for j = 1:nPerGroup
                    offset = (j-(nPerGroup+1)/2)*(boxHeight*1.3);
                    rectangle(axBottom, 'Position',[Q1_low(i,j), yBase+offset-boxHeight/2, Q3_low(i,j)-Q1_low(i,j), boxHeight],...
                              'FaceColor',greyRGBs(j,:),'EdgeColor','none');
                    % Median dot
                    if add_color_legend
                        plot(axBottom, MED_low(i,j), yBase+offset, 'o', ...
                        'MarkerSize',5, ...
                        'MarkerFaceColor', cmap(scaledIdx_low(i,j),:), ...
                        'MarkerEdgeColor', cmap(scaledIdx_low(i,j),:));
                    else
                        plot(axBottom, MED_low(i,j), yBase+offset,'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
                    end
                end
            end

        
            yticks(axBottom, 1:nGroups*groupSep)
            yticklabels(axBottom, strrep(rxn_names_low,"_","\_"))
            xlabel(axBottom,'Flux value [mMol/(gDW*h)]')
        end

        %%% =========================
        % Set axes limits based on MED only
        %%% =========================
        if exist('axTop','var')
            % X limits
            xMinHigh = min(MED_high(:));
            xMaxHigh = max(MED_high(:));
            % Add small padding
            xPad = (xMaxHigh - xMinHigh) * 0.05;
            %xlim(axTop, [xMinHigh-xPad, xMaxHigh+xPad]);
        
            % Y limits
            yMinHigh = 0.5;  % first y
            yMaxHigh = nHigh * groupSep + 0.5;
            ylim(axTop, [yMinHigh, yMaxHigh]);
        end
        
        if exist('axBottom','var')
            % X limits
            xMinLow = min(MED_low(:));
            xMaxLow = max(MED_low(:));
            xPad = (xMaxLow - xMinLow) * 0.05;
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
        if add_color_legend
            % Choose axTop if it exists, else axBottom
            if exist('axTop','var')
                cbAx = axTop;
            else
                cbAx = axBottom;
            end
            
            colormap(cbAx, cmap);          % Set colormap for chosen axis
            % --- Set caxis safely
            if valMin == valMax
                caxis(cbAx, [valMin-0.5, valMax+0.5]);  % small padding if single value
            else
                caxis(cbAx, [valMin, valMax]);
            end

            
            cb = colorbar(cbAx);           % Attach colorbar
            cb.Label.String = ...
               'Reduced Cost';
            cb.FontSize = 14;
        end
        
        % --- Grey patches legend (for models)
        % Choose axTop if it exists, else axBottom
        if exist('axTop','var')
            lgdAx = axTop;
        else
            lgdAx = axBottom;
        end
        
        hGrey = gobjects(nPerGroup, 1);
        for j = 1:nPerGroup
            hGrey(j) = patch(lgdAx, NaN, NaN, greyRGBs(j,:), 'EdgeColor', 'none');
        end
        
        % Add legend
        lgd = legend(lgdAx, hGrey, modelList, 'Location','northeastoutside');
        % lgd = legend(lgdAx, hGrey, modelList, 'Location','northeast');
        lgd.FontSize = 16;
        lgd.Box = 'off';
        lgd.Color = 'none';
        lgd.Title.String = "Models";

    end

    % --- Reverse direction if needed ---
    if options.thresholdFlux=="upper"
        set([axTop axBottom],'XDir','reverse')
    end
    

    %%% =========================
    % Table
    % =========================
    % resort the table according to the order in the plots
    if exist('T', 'var')
        T = T([flip(string(rxn_names_high));flip(string(rxn_names_low))],:);
        tbl = uitable(fig, ...
                        'Data', T, ...
                        'ColumnName', T.Properties.VariableNames, ...
                        'Units','normalized', ...
                        'Position',[plotWidth+0.05 0.10 0.40 0.85], ... % width increased from 0.30 → 0.40
                        'FontSize',16, ...
                        'ColumnWidth','auto');
    
        
         tbl.FontSize = 16; 
    end

end




