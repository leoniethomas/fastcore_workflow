function fig = get_flux_plot(project,comparison_name, idx_to_vis,options)
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
    % - idx_to_vis:             positions of the rxns to be displayed in the
    %                           choosen reference model
    % - options:                - FVA = true (default false) to display the FVA boundaries
    %                             around the FBA solution as a grey box.
    %                           - reducedCost = true (default false) to display the reduced
    %                             cost values for every reaction as the color
    %                             of the FBA dot.
    %                           - threshold_flux= wether to apply an upper
    %                             lower or no threshold to the selected
    %                             reaction fba values
    %                           - 
    % 
    % Output:                   Display of figure                
    % 

    arguments
        project
        comparison_name (1,1) string
        idx_to_vis
        options.FVA  (1,1) logical = false
        options.reducedCost (1,1) logical = false
        options.threshold_flux (1,1) string {mustBeMember(options.threshold_flux,["lower", "upper","none"])} ="none" 
        options.title_plots = ""
    end


    modelList = project.comparisons.(comparison_name).modelNames;
    reference_model = project.comparisons.(comparison_name).reference_model;
    
    replacement_value = "analysis.FBA.v"; % get the fba solution values
    ordered_fba_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
    
    
    if options.threshold_flux == "upper"
        get_exchange_rxns_idx = intersect(find(sum(ordered_fba_matrix,2) ~=0 & mean(ordered_fba_matrix,2) <0), ...
                                          idx_to_vis);  
        title_word = "positive flux reactions";
    elseif options.threshold_flux == "lower"
        get_exchange_rxns_idx = intersect(find(sum(ordered_fba_matrix,2) ~=0 & mean(ordered_fba_matrix,2) > 0), ...
                                                      idx_to_vis); 
        title_word = "negative flux reactions";
    elseif options.threshold_flux == "none"
        get_exchange_rxns_idx = intersect(find(sum(ordered_fba_matrix,2) ~=0), ...
                                                      idx_to_vis);
        title_word = options.title_plots;
    else
        error("wrong value for threshold choosen. Possible values: lower, upper, none")
    end
    ordered_fba_matrix_ex = ordered_fba_matrix(get_exchange_rxns_idx,:);

    rxn_names = project.models.(reference_model).model.rxns(get_exchange_rxns_idx);
    rxn_formulas = printRxnFormula(project.models.(reference_model).model, 'rxnAbbrList', rxn_names, 'printFlag', false);

    % medium composition 
    % do the models have the same  medium ? 
    % if so I can add it as column to the plot, otherwise it is not added
    models = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    media_for_models = structfun(@(x) x.settings.medium, models);
    ref = fieldnames(models);
    ref = models.(ref{1,1}).settings.medium;
    medium_is_equal_between_models = all(arrayfun(@(x) isequaln(x, ref), media_for_models));
    

    if medium_is_equal_between_models
        medium_constrained = ismember(rxn_names, ref.medium_composition.ExRxns_Recon3D) | ...
                         ismember(rxn_names, ref.manual_set_boundaries.unwanted_export) | ...
                         ismember(rxn_names, ref.manual_set_boundaries.unwanted_import);
        

        ordered_lb = getOrderedFeatureMatrix(project,"consistent_medium_constrained_model",...
                                             "rxns",reference_model,"model.lb");
        ordered_ub = getOrderedFeatureMatrix(project,"consistent_medium_constrained_model",...
                                             "rxns",reference_model,"model.ub");
        ordered_ub = ordered_ub(get_exchange_rxns_idx,:);
        ordered_lb = ordered_lb(get_exchange_rxns_idx,:);

        T = table(rxn_formulas, medium_constrained, ...
                  'VariableNames', {'Reaction Formula','medium constrained'}, ...
                  'RowNames',rxn_names);
        T = T(flip(string(T.Properties.RowNames)),:);
        T.lb = flip(ordered_lb);
        T.ub = flip(ordered_ub);
    end
    


    if options.FVA
        replacement_value = "analysis.FVA.maxFlux"; % get the fba solution values
        ordered_fva_max_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
        replacement_value = "analysis.FVA.minFlux"; % get the fba solution values
        ordered_fva_min_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
        ordered_fva_max_matrix_ex = ordered_fva_max_matrix(get_exchange_rxns_idx,:);
        ordered_fva_min_matrix_ex = ordered_fva_min_matrix(get_exchange_rxns_idx,:);
    end
    if options.reducedCost
        replacement_value = "analysis.FBA.basis.reducedcost"; % get the fba solution values
        ordered_reducedCost_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
        ordered_reducedCost_matrix_ex = ordered_reducedCost_matrix(get_exchange_rxns_idx,:);
    end
    %if options.shadowPrice
    %    replacement_value = "analysis.FBA.basis.dual"; % get the fba solution values
    %    % shadow prices are measured for every metabolite therefore mapped according to the mets field
    %    ordered_shadowPrices_matrix = getOrderedFeatureMatrix(project,modelList,"mets",reference_model,replacement_value);
    %end
    
    if ~options.FVA & ~options.reducedCost % specify which of the fields need to be true and false!!!
        % in the case that only the FBA solution should be visualized we
        % use a grouped horizontal barplot to do so
       
        nRxns = size(ordered_fba_matrix_ex,1);
        
        % Create a UI figure
        fig = uifigure('Name', title_word + " with Table",'Position',[100 100 1000 450]);
        
        % Left plot area width (normalized)
        plotWidth = 0.65;  % 65% for plot, 35% for table
        
        % Create axes for horizontal grouped bar plot
        ax = uiaxes(fig);
        %ax.Layout.Tile = [];  % not using grid, manual positioning below
        ax.Units = 'normalized';
        ax.Position = [0.02 0.1 plotWidth 0.85];  % left, bottom, width, height
        
        % Plot horizontal grouped bars
        barh(ax, ordered_fba_matrix_ex, 'grouped');  
        grid(ax,'on')
        ax.GridColor = [0.8 0.8 0.8];
        ax.GridAlpha = 0.5;
        if options.threshold_flux == "upper"
            ax.XDir = 'reverse';
        end
        ax.FontSize = 14;
        
        % Labels and title
        yticks(ax, 1:nRxns)
        yticklabels(ax, strrep(rxn_names, "_", "\_"))
        xlabel(ax, 'Flux value [mMol/(gDW*h)]')
        title(ax, title_word)
        legend(ax, modelList, 'Location','best')
        
        % --- Create UITable beside the bar chart ---
        tbl = uitable(fig, ...
            'Data', T, ...
            'ColumnName', T.Properties.VariableNames, ...
            'Units','normalized', ...
            'Position',[plotWidth+0.05 0.1 0.30 0.85]);  % position to the right
        
        tbl.FontSize = 14;


    else
        Q1  = ordered_fva_min_matrix_ex; %FVAmin
        MED = ordered_fba_matrix_ex; %FBA sol
        Q3  = ordered_fva_max_matrix_ex; %FVAmax
    
        if options.reducedCost
            add_color_legend = 1;
            cmap = colormap('cool')
            N = size(cmap,1);
            % check if you have any values for the reduced cost 
            if length(unique(ordered_reducedCost_matrix_ex)) == 1
                % if there is only one value in the dataframe just set
                % the scaledIdx to be that value
                scaledIdx  = ordered_reducedCost_matrix_ex + 1; 
                valMin = unique(ordered_reducedCost_matrix_ex);
                valMax = unique(ordered_reducedCost_matrix_ex) + 1;
                
            else
                % Example variable to color medians (normalized 0–1)
                
                valMin = min(ordered_reducedCost_matrix_ex(:));   % e.g., 0
                valMax = max(ordered_reducedCost_matrix_ex(:));   % e.g., ~11.45 in your data
                % Map each value into [1, N] range
                scaledIdx = round( (ordered_reducedCost_matrix_ex - valMin) / (valMax - valMin) * (N-1) ) + 1;
                
                % Clip to valid range
                scaledIdx(scaledIdx < 1)   = 1;
                scaledIdx(scaledIdx > N)   = N;
                
            end

        else
            add_color_legend = 0;
        end

        %%

        
        % --- Create UI figure and UIAxes
        fig = uifigure('Name','Plot Left + Table Right','Position',[100 100 900 400]);
        % Size settings
        margin     = 20;     % space around components
        tableWidth = 300;    % width allocated for the table
        figWidth   = fig.Position(3);
        figHeight  = fig.Position(4);
        
        plotWidth  = figWidth - tableWidth - 3*margin;  % remaining width
        plotHeight = figHeight - 2*margin;              % full height minus margins

        
        % Set up axes on left (adjust sizes as needed)
        ax = uiaxes(fig, ...
                    'Position', [20 50 550 330]);
        
        % Important: hold on *the UIAxes*, not default axes
        hold(ax, 'on')
        
        % Plot your custom graphic on the UIAxes
        boxHeight = 0.2;
        groupSep  = 1;
        [nGroups,nPerGroup] = size(Q1);
        
        greyLevels = linspace(0.9, 0.6, nPerGroup);
        greyRGBs = [greyLevels', greyLevels', greyLevels'];
        
        for i = 1:nGroups
            yBase = i * groupSep;  % base y for this category
            for j = 1:nPerGroup
                offset = (j - (nPerGroup+1)/2) * (boxHeight * 1.3);
        
                % Draw grey box
                rectangle(ax, ...
                  'Position',[Q1(i,j), yBase+offset - boxHeight/2, ...
                              Q3(i,j)-Q1(i,j), boxHeight], ...
                  'FaceColor', greyRGBs(j,:), ...
                  'EdgeColor', 'none');
                
                if add_color_legend
                    % Plot median dot
                    plot(ax, MED(i,j), yBase+offset, 'o', ...
                         'MarkerSize', 5 , ...
                         'MarkerFaceColor', cmap(scaledIdx(i,j),:), ...
                         'MarkerEdgeColor', cmap(scaledIdx(i,j),:));
                else
                    plot(ax, MED(i,j), yBase+offset, 'o', ...
                         'MarkerSize', 5 , ...
                         'MarkerFaceColor', 'k', ...
                         'MarkerEdgeColor', 'k');
                end
            end
        end

        % --- Set ticks and labels on UIAxes
        yticks(ax, (1:nGroups)*groupSep)
        yticklabels(ax, strrep(rxn_names, "_", "\_"))
        
        xlabel(ax, 'Flux Value [mMol/(gDW*h)]', 'FontSize', 16)
        title(ax, "FVA, FBA + Reduced cost for the " + title_word, 'FontSize', 18)
        
       
        
        % Reverse x direction (if intended)
        if options.threshold_flux == "upper"
            ax.XDir = 'reverse';
            xlim(ax, [-8 0.1])
        else
            xlim(ax, [-0.1 8])
        end
        
        % Apply grid
        grid(ax, 'on')
        ax.GridColor = [0.8 0.8 0.8]
        ax.GridAlpha = 0.5
        ax.FontSize = 16
        
        % --- Add colorbar if needed
        if add_color_legend
            colormap(ax, cmap)            % colormap for UIAxes
            caxis(ax, [valMin valMax])    % set real data range for color mapping
            cb = colorbar(ax)             % attach to UIAxes
            cb.Label.String = ...
               'Sensitivity of objective function to changes in the flux boundaries'
            cb.FontSize = 14;
        end
        
        
        hGrey = gobjects(nPerGroup, 1);
        for j = 1:nPerGroup
            hGrey(j) = patch(ax, NaN, NaN, greyRGBs(j,:), 'EdgeColor', 'none');
        end
        
        % Add legend *on UIAxes*
        lgd = legend(ax, hGrey, modelList, 'Location', 'northeastoutside');
        lgd.FontSize = 16;
        lgd.Title.String = "Models";

        %%
        % Create uitable ON THE RIGHT
        tbl = uitable(fig, ...
                      'Data', T, ...
                      'ColumnName', T.Properties.VariableNames, ...
                      'Position', [plotWidth + 2*margin, margin, tableWidth, plotHeight]);
        
        % Optional: adjust font sizes
        tbl.FontSize = 14;
        ax.FontSize = 12;
    
    end
    

end

