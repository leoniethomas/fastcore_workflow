function project = modelsComparison(project, modelList,reference_model,analyses,identifier)
    % This function runs a set of analysis for the comparison of the
    % specified genes.
    % A number of analysis are run: 
    % - structural analysis: based on the differential presence of
    %   metabolites, genes and reactions in the different models 
    % - functional analysis: based on the quantitative values like FVA,
    %   FBA, sampling 
    % 
    % Inputs: 
    %   - project:          the object which is the output of the single_model_analysis
    %                       entailing the results of fba,fva,sampling, single gene
    %                       deletion etc. for a single model 
    %   - modelList:        the list of Model names to be included in the comparison 
    %   - reference_model:  the reference model used to compute the relative feature presence to 
    %   - analyses:         the list of analyses which should be performed 
    %   - identifier:       a string, will be added as a postfix to the analysis name, can be choosen freely
    %                       default)
    %   
    % Output : 
    %   - project:          project object with a added comparison field entailing
    %                       all the output, modelcomparison information

    arguments
        project
        modelList (1,:) string = strings(0)
        reference_model (1,:) string = "orig_model"
        analyses  (1,:) string = ["modelStructureComparison"]
        identifier (1,:) string = string(datetime('now','Format','_yyyyMMdd_HHmmss'))
    end
    
    % if no list of genes is given just compare all the models in the
    % project.models slot
    if isempty(modelList)
        modelList = project.models;
    end
    
    % give the comparison the name of all models + a identifier choosen
    comparison_name = join(modelList, "_vs_") + identifier;

    % run structural model comparison
    project.comparisons.(comparison_name) = modelStructuralComparison(project,modelList,reference_model);
    project.comparisons.(comparison_name).reference_model = reference_model; 
    % run functional model comparison


    modelFunctionalComparison(project, comparison_name,analyses);
    

end

function modelFunctionalComparison(project, comparison_name,analyses)
    % The structure comparison is a function that compares the models
    % listed on structural differences. Structural differences in the
    % context of Fastcore can be defined as the set of reactions that are
    % kept when running fastcore. This means we check for the existence of
    % rxns, metabolites and genes in the model, and their overlap between
    % models.

    arguments
        project
        comparison_name
        analyses
    end

    modelList = project.comparisons.(comparison_name).modelNames;
    reference_model = project.comparisons.(comparison_name).reference_model;
    % first visualize the fba solution 
    % what are the variables I need from the fba solution 
    
    fba_objective_values = cell2mat(cellfun(@(x) project.models.(x).analysis.FBA.f(1,1) ,modelList,"UniformOutput",false));
    get_exchange_rxns_idx = find(findExcRxns(project.models.(reference_model).model));    
    
    figure
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

    %%% barplot for the biomass/defined objective 

    nexttile(1);
    bar(fba_objective_values)
    title('Model comparison: flux of optimized reaction')
    ylabel('Reaction flux value for objective function [mMol/(gDW*h)]')
    xlabel('Model')
    xticklabels(modelList)
    
    get_flux_plot(project, comparison_name,get_exchange_rxns_idx)

    % compute Fluxvariability display in heatmap 
    %%% -> model wide 

    %%% -> per subsystem

    % compute Fluxsum 

    %%% -> show the top 20 most variant metabolites excluding known cofactors 
    % cofactorNames = ["atp", "adp", "amp", "nad", "nadh", "nadp", "nadph", ...
    %             "coa", "accoa", "fad", "fadh2", "pi", "pp_i"];

    %%% -> show the fluxsum in the defined pathways
    % get_essential_pathway_metabolites('Glycolysis',project,reference_model)

end


function structure_analysis = modelStructuralComparison(project, modelList,reference_model)
    % The structure comparison is a function that compares the models
    % listed on structural differences. Structural differences in the
    % context of Fastcore can be defined as the set of reactions that are
    % kept when running fastcore. This means we check for the existence of
    % rxns, metabolites and genes in the model, and their overlap between
    % models.

    arguments
        project
        modelList
        reference_model
    end

    % extract models you want to compare from the project object
    models_list = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    models = structfun(@(x) x.model, models_list, 'UniformOutput',false);
    structure_analysis.modelNames = string(fieldnames(models));


    % get the size of the different models
    data = struct2array(structfun(@(x) {numel(x.rxns); numel(x.mets); numel(x.genes)}, ...
                           models, 'UniformOutput', false))';
    array2table(data, ...
                    'VariableNames', {'count_reactions','count_metabolites','count_genes'}, ...
                    'RowNames', string(fieldnames(models))')
   
    %%%% get the indices of all the rxns from the different models
    %%%% specified to put them into one table

    
    [~,rxn_mapping] = getOrderedFeatureMatrix(project,modelList,"rxns", reference_model);
    structure_analysis.rxn_mapping_table = array2table(rxn_mapping,"VariableNames",modelList,"RowNames",string(project.models.(reference_model).model.rxns))
    
    %%%%%%%%%%%%%%%%%%%% get the intersections/outersection 
    %%%%%% similarities between models based on rxns,genes, mets

    for field_to_investigate = ["genes", "mets", "rxns"]
        [ordered_feature, ~] = getOrderedFeatureMatrix( ...
            project, modelList, field_to_investigate, reference_model);
    
        % Plot Venn / Heatmap of intersections based on presence
        plotFlexibleVenn( ...
            ordered_feature, ...
            structure_analysis.modelNames, ...
            "Structural model comparison: " + field_to_investigate + " presence");
        % plot Venn/Heatmap of the overall intersection/outersection of the
        % models based on rxn prescence 
        plotFlexibleVenn(ordered_feature,...
                             structure_analysis.modelNames, ... 
                             "Structural model comparison: " + field_to_investigate + " presence");
    
        % get the jaccard distances - based on reaction presence
        Jacc_distance = 1 - squareform(pdist(ordered_feature','jaccard'));
        title_fig = "Jaccard similarity of " + field_to_investigate + " presence (0 or 1) between models ";
        figure
        heatmap(structure_analysis.modelNames,structure_analysis.modelNames, Jacc_distance);
        title(title_fig);
    end
    
    %%%%%%%%%%%%%%%%%%%%% get presence of rxns in each subsytem/pathway
        
    pathways = string(project.models.(reference_model).model.subSystems); % get pathways from reference model
    unique_pathways = unique(pathways); 

    % for every pathway get the matrix identifying the rnxs from reference
    % model in this pathway
    groups = arrayfun(@(x) find(pathways == x), unique_pathways, 'UniformOutput', false);
    num_groups = length(groups);
    G = zeros(num_groups, size(ordered_feature,1));

    for g = 1:num_groups
        G(g, groups{g}) = 1;
    end
    
    % get the presence per subsystem in the context specific models 
    % by mulitplying the rxns presence for each subsystem (matrix ordered features) 
    % with the matrices defining the subsystem for every rxns
    pathway_counts = array2table(G * ordered_feature, ...
                                 'VariableNames', structure_analysis.modelNames,...
                                 'RowNames',cellstr(unique_pathways));
    % add reference model count to be able to make a relative abundance
    pathway_counts.reference_model = groupcounts(pathways);
    
    relative_counts = array2table(pathway_counts{:,1:end-1} ./ pathway_counts.reference_model, ...
                                 'VariableNames', structure_analysis.modelNames,...
                                 'RowNames',cellstr(unique_pathways));
    
    % get the idx of the most variant pathways in terms of rxns presence
    relative_counts.row_var = var(relative_counts{:,:}, 0, 2);
    pathway_counts.row_var = var(pathway_counts{:,1:end-1},0,2);
    pathway_counts = pathway_counts(pathway_counts.reference_model < 1000,:);
    % Get indices of top n highest variance rows

    
    pathway_counts = sortrows(pathway_counts,"row_var","descend");
    pathway_counts = pathway_counts(find(pathway_counts.row_var ~= 0),:);
    relative_counts = relative_counts(pathway_counts.Properties.RowNames,:);
    % plot top 20 most variant pathways between the choosen models
    data = relative_counts{:,1:end-1};
    rowNames = string(relative_counts.Properties.RowNames);
    colNames = structure_analysis.modelNames;

    figure
    tiledlayout(1,4, ...
        'TileSpacing','compact', ...
        'Padding','compact')
    
    % ---- Bar plot (LEFT) ----
    ax1 = nexttile(1);
    barh(pathway_counts.reference_model)
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
    ax2 = nexttile(2,[1 3]);
    
    imagesc(data)
    nColors = 256;
    whiteToBlue = [linspace(1,0,nColors)', linspace(1,0,nColors)', ones(nColors,1)];
    colormap(ax2, whiteToBlue)
    colorbar
    title("relative counts of subsystem rxn occurence/reference model" )    % grayscale
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
            % (or multiply relative_counts by reference_model if needed)
            value = pathway_counts{i,j}; % +1 because first column is reference_model
            % Place text at the center of the tile
            text(ax2, j, i, num2str(value), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'Color','k', ...          % black text
                'FontSize',10)
        end
    end

    
    % ---- Align rows ----
    linkaxes([ax1 ax2],'y')

    %%% Visualize the core reaction per model
    data = struct2cell(structfun(@(x) [ length(x.core_reactions) - sum(ismember(x.core_reactions, x.model.rxns)); ...
                                            sum(ismember(x.core_reactions, x.model.rxns));...
                                            length(x.model.rxns) - sum(ismember(x.core_reactions, x.model.rxns));...
                                            length(x.model.rxns)], ...
                                            models_list, 'UniformOutput', false));
    data = [data{:}];
    
    % ---- Create layout ----
    upper_data = data(2:3,:);

    %figure
    categories = fieldnames(models_list)';  % model names
    %t = tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact');

    figure
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
    
    % --- first barplot
    ax1 = nexttile(1);
    bar(upper_data', 'stacked')
    
    
    % Labels
    set(gca, 'XTickLabel', categories, 'FontSize', 14)  % increase tick labels font
    xlabel('Models', 'FontSize', 14)
    ylabel('# rxns', 'FontSize', 14)
    
    % Legend
    legend({"non-core reactions", "core reactions"}, 'Location','northwest', 'FontSize', 14)
    
    % Title
    title('Core and non-core reactions per model', 'FontSize', 14)

    % ---- Compute percentages of upper stack ----
    total = sum(upper_data,1);                   % total per model
    percent_upper = 100 * upper_data(2,:) ./ total;  % percentage of upper bar
    
    % ---- Add text labels on top of upper bars ----
    for i = 1:size(upper_data,2)  % loop over models
        % x-position is the bar center, y-position is height of lower + upper
        x = i;
        y = upper_data(1,i) + upper_data(2,i)/2;  % middle of upper stack
        text(x, y, sprintf('%.1f%%', percent_upper(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold');
    end

    
    
    
     % --- second barplot
     data = data(1:2,:);  % only core vs non-core counts
    ax2 = nexttile(2);
    hb = bar(data', 'stacked');
    
    % Labels
    set(gca, 'XTickLabel', categories, 'FontSize', 14)
    xlabel('Models', 'FontSize', 14)
    ylabel('# rxns', 'FontSize', 14)
    
    % Legend
    legend({ "not included","included"}, 'Location','northwest', 'FontSize', 14)
    
    % Title
    title('Core reactions used to construct the models.', 'FontSize', 14)
    
    % ---- Compute percentages of upper stack ----
    total = sum(data,1);                   % total per model
    percent_upper = 100 * data(2,:) ./ total;  % percentage of upper bar
    
    % ---- Add text labels on top of upper bars ----
    for i = 1:size(data,2)  % loop over models
        % x-position is the bar center, y-position is height of lower + upper
        x = i;
        y = data(1,i) + data(2,i)/2;  % middle of upper stack
        text(x, y, sprintf('%.1f%%', percent_upper(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold');
    end

    % Tile 3 (THIS IS THE KEY)
    ax3 = nexttile(3,[1 2]);   % column 2, span both rows
    axis(ax3,'off')
    hold(ax3,'on')

    core_reactions_included = struct2cell(structfun(@(x) x.core_reactions(find(ismember(x.core_reactions, x.model.rxns)))', ...
                                            models_list, 'UniformOutput', false));
    core_reactions_included = unique([core_reactions_included{:}]);
    
    core_presence = structure_analysis.rxn_mapping_table{core_reactions_included,:} ~= 0;
    figV = plotFlexibleVenn(core_presence, structure_analysis.modelNames, ... 
                     "Structural model comparison: core rxns presence");

    % Find axes inside Venn figure
    axV = findobj(figV,'Type','axes','-not','Tag','legend');
    axV = axV(1);
    
    % Copy graphics
    copyobj(allchild(axV), ax3)

    
    
    % Fix geometry
    axis(ax3,'tight')
    axis(ax3,'equal')
    ax3.Clipping = 'off';

    
    close(figV)


end


function [ordered_feature_matrix,ordered_rxn_matrix_idx] = getOrderedFeatureMatrix(project,modelList,field_to_investigate,reference_model,replacement_value)
    % this function brings the features for multiple models (for example the rxns presence (0
    % or 1)) into the same order than the reference model specified. 
    % Input: 
    % - project:                project object which is the output of
    %                           single_model_analysis script
    % - modelList:              names from models from which to get the
    %                           feature precence
    % - field_to_investigate:   feature presence of interest as string
    %                           "genes", "mets", or "rxns"
    % - reference_model:        reference model that give the order after
    %                           which the choosen feature will be ordered 
    % - replacement_value:      which value to put into the matrix. Just a
    %                           1 indicating that the given feature is present 
    %                           or not (0) in each of the choosen models, 
    %                           or the actuall values of sampling,fba, or fva ?  
    % 
    % Output: 
    % - ordered_feature_matrix: Matrix storing the wanted features 
    %                           (presence of genes,mets,rxns,or actuall 
    %                           fva,fba,sampling values) in the order of 
    %                           the reference model.
    %                           dim nxm with n = features in reference model
    %                                   and  m = models to compare
    % 
    % - ordered_rxn_matrix_idx: Tabel with the actual positions the rxns are 
    %                           stored at in the individual models. 
    %                           dim nxm with n = rxns in reference model
    %                                   and  m = models to compare
    % 
    arguments
        project
        modelList
        field_to_investigate (1,1) string ="rxns"
        reference_model  = strings(0)
        replacement_value (1,1) =1
    end
    % by default we bring all the models in the same order as the original model
    if isempty(reference_model)
        reference_model = project.models.orig_model.model;
    else
        reference_model = project.models.(reference_model).model;
    end

    models = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    %models = structfun(@(x) x.model, models, 'UniformOutput',false);

    ordered_feature_matrix = struct2array(structfun(@(x) getOrderedFeature(x,reference_model,field_to_investigate,replacement_value), ...
                                             models,'UniformOutput',false));
    
    ordered_rxn_matrix_idx = struct2array(structfun(@(x) getOrderedFeature(x,reference_model,"rxns","idx"), ...
                                             models,'UniformOutput',false));
end

function ordered_feature = getOrderedFeature(model,reference_model,field_to_investigate,replacement_value)
    % this gives you back a vector in the length of the defined reference
    % data field (for example the length of rxns in the reference model)
    % and then give a 1 if this feature (for example the rxns) is present
    % in the specified model and a 0 if it is not present in the specified
    % model.
    % this function brings the features for a single model(for example the rxns presence (0
    % or 1)) into the same order than the reference model specified. 
    % Input: 
    % - model:                  model from which to obtain the feature
    % - reference_model:        reference model that give the order after
    %                           which the choosen feature will be ordered 
    % - field_to_investigate:   feature presence of interest as string
    %                           "genes", "mets", or "rxns"
    % - replacement_value:      which value to put into the vector. Just a
    %                           1 indicating that the given feature is present 
    %                           or not (0) in each of the choosen models, 
    %                           or the actuall values of sampling,fba, or fva ?  
    % 
    % Output: 
    % - ordered_feature:        a vector storing the wanted features 
    %                           (presence of genes,mets,rxns,or actuall 
    %                           fva,fba,sampling values) 
    %                           in the order of the reference model.
    %                           length(ordered_feature) == n with 
    %                           n = length(model.(field_to_investigate))
    % 
    arguments
        model
        reference_model
        field_to_investigate (1,1) string ="rxns"
        replacement_value (1,1) =1
    end

    ordered_feature_idx = zeros(size(reference_model.(field_to_investigate),1),size(reference_model.(field_to_investigate),2));
    [~,idx] = ismember(model.model.(field_to_investigate),reference_model.(field_to_investigate));
    ordered_feature_idx(idx) = 1:numel(model.model.(field_to_investigate));
    
    % now we have a matrix that stores the idx of the features of the
    % models on the position of where the feature is in the reference
    % model.
    ordered_feature = ordered_feature_idx;
    if isnumeric(replacement_value)
        ordered_feature(ordered_feature_idx >0) = replacement_value;
    elseif isstring(replacement_value) || ischar(replacement_value)
        if replacement_value ~= "idx"
            replacement_value = strsplit(replacement_value,".");
            replacement_values = model.(replacement_value(1));
            for x = replacement_value(2:end)
                replacement_values = replacement_values.(x);
            end
            if length(replacement_values) ~= length(model.model.(field_to_investigate))
                error("The size of the slot you have choosen (rxns, genes) does not aggree with the replacement value slot! Are you sure the replacement slot you choose really contains information about the choosen slot (genes, rxns,mets) ??")
            end
            ordered_feature(idx) = replacement_values;
        end
    end

end


function fig = get_flux_plot(project,comparison_name, idx_to_vis,options)

    arguments
        project
        comparison_name
        idx_to_vis
        options.FVA  (1,1) logical = true
        options.shadowPrice (1,1) logical = false
        options.reducedCost (1,1) logical = true
        options.threshold_flux = ["import", "export"]
        options.title_plots = ""
    end


    modelList = project.comparisons.(comparison_name).modelNames;
    reference_model = project.comparisons.(comparison_name).reference_model;
    
    replacement_value = "analysis.FBA.v"; % get the fba solution values
    ordered_fba_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
    
    for threshold = options.threshold_flux
        if threshold == "import"
            get_exchange_rxns_idx = intersect(find(sum(ordered_fba_matrix,2) ~=0 & mean(ordered_fba_matrix,2) <0), ...
                                              idx_to_vis);  
            title_word = "import reactions";
        elseif threshold == "export"
            get_exchange_rxns_idx = intersect(find(sum(ordered_fba_matrix,2) ~=0 & mean(ordered_fba_matrix,2) > 0), ...
                                                          idx_to_vis); 
            title_word = "export reactions";
        elseif threshold == "all"
            get_exchange_rxns_idx = intersect(find(sum(ordered_fba_matrix,2) ~=0), ...
                                                          idx_to_vis);
            title_word = options.title_plots;
        else
            error("wrong value for threshold choosen. Possible values: import or export")
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
            %T = T(flip(project.models.(reference_model).model.rxns(get_exchange_rxns_idx)),:);
            T = T(flip(string(T.Properties.RowNames)),:);
            T.lb = flip(ordered_lb);
            T.ub = flip(ordered_ub);
            
        else
            
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
            if threshold == "import"
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
                    cmap = colormap('cool')
                    N = size(cmap,1);
                    
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
            if threshold == "import"
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

end



function essential_pathway_metabolites = get_essential_pathway_metabolites(pathway,project,reference_model)
    % source for pathway metabolite composition: https://github.com/sysbiolux/MidbrainOrganoid_Miro1_scMetabMod/blob/main/ModelAnalysis/compareFBA_v3_PAPER.m
    
    arguments
        pathway (1,1) string {mustBeMember(pathway,["Glycolysis","PPP","TCA_cytoplasm_only","Tricarboxylic_acid_cycle","Oxidative_phosphorylation","Fatty_acid_oxidation","Cholesterol","Dopamine_Tyrosine_metabolism","Fatty_acid_oxidation__all__crn","Fatty_acid_oxidation__all__coa"])}
        project
        reference_model
    end
        pathways.Glycolysis = ["glc_D[c]","g6p[c]","f6p[c]","fdp[c]","dhap[c]","g3p[c]","13dpg[c]","3pg[c]","2pg[c]","pep[c]","pyr[c]"];
        pathways.PPP=["ru5p_D[c]","s7p[c]","e4p[c]"];
        pathways.TCA_cytoplasm_only= ["cit[c]","icit[c]","akg[c]","succ[c]","fum[c]","oaa[c]"];
        pathways.Tricarboxylic_acid_cycle=["cit[m]","icit[m]","akg[m]","succoa[m]","succ[m]","fum[m]","mal_L[m]","oaa[m]"];
        pathways.Oxidative_phosphorylation=["nadh[m]","fadh2[m]","focytC[m]","q10h2[m]","atp[m]"];
        pathways.Fatty_acid_oxidation=["accoa[c]","accoa[m]","coa[c]","coa[m]"];
        %%
        temp=find(ismember(string(project.models.(reference_model).model.subSystems),'Fatty acid oxidation'));
        metList=findMetsFromRxns(project.models.(reference_model).model, ...
                                 project.models.(reference_model).model.rxns(temp))';
        temp=find(contains(metList,'coa'));
        pathways.Fatty_acid_oxidation__all__coa=string(metList(temp));
        %%
        temp=find(ismember(string(project.models.(reference_model).model.subSystems),'Fatty acid oxidation'));
        metList=findMetsFromRxns(project.models.(reference_model).model, ...
                                 project.models.(reference_model).model.rxns(temp))';
        temp=find(contains(metList,'crn'));
        pathways.Fatty_acid_oxidation__all__crn=string(metList(temp));
        %%
        pathways.Cholesterol=["sql[r]","zymst[r]","chsterol[r]","chsterol[c]"];
        pathways.Dopamine_Tyrosine_metabolism=["phe_L[c]","tyr_L[c]","34dhphe[c]","dopa[c]","34dhpac[c]","34dhpha[c]","34dhpe[c]","tym[c]"];

        essential_pathway_metabolites = pathways.(pathway);
end




