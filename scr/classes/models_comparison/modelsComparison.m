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
    project.comparisons.(comparison_name) = modelFunctionalComparison(project, comparison_name,analyses);


end

function structure_analysis = modelFunctionalComparison(project, comparison_name,analyses)
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
    field_to_investigate = "rxns"; % the fba solution values are measurements belonging to the rxns
    reference_model = project.comparisons.(comparison_name).reference_model;
    % first visualize the fba solution 
    % what are the variables I need from the fba solution 
    
    replacement_value = "FBA.v"; % get the fba solution values
    ordered_fba_matrix = getOrderedFeatureMatrix(project,modelList,field_to_investigate,reference_model,replacement_value);
    replacement_value = "FBA.basis.reducedcost"; % get the fba solution values
    ordered_reducedCost_matrix = getOrderedFeatureMatrix(project,modelList,field_to_investigate,reference_model,replacement_value);
    replacement_value = "FBA.basis.dual"; % get the fba solution values
    field_to_investigate = "mets"; % shadow prices are measured for every metabolite therefore mapped according to the mets field
    ordered_shadowPrices_matrix = getOrderedFeatureMatrix(project,modelList,field_to_investigate,reference_model,replacement_value);



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

    ordered_feature_idx = zeros(length(reference_model.(field_to_investigate)),1);
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
            replacement_values = model.analysis.(replacement_value(1));
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




