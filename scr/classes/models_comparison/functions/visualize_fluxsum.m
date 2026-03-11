function fluxsum_sets = visualize_fluxsum(project,comparison_name,met_idx,rxn_idx,rxn_set_labels,plot_type,exclude_coenzymes,ignore_compartment)
    % This function visualizes the fluxsum either in a heatmap or in a
    % violinplot between different models to enable easy comparison of
    % model results in sampling. 
    % Before running this script a modelsComparison must be run. Then the
    % name of the comparison performed (given as output of the
    % modelsComparison function) can be given as an input together with the
    % project object entailing the comparison + all the needed models. 
    %
    % The function produces as many figures as you define sets in the
    % rxn_idx argument. This is meant to enable you to visualize different
    % sets/subsystems of your choosing to be visualized in one Figure. 
    % 
    % The function allows visualization of the fluxsum on different levels.
    % If you choose the violin plot_type -> the fluxsum for every metabolite is
    % computed and then visualized in a separate violinplot for every
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
    %       - comparison_name: name of the comparison you want to visualize
    %         (the second output of the modelsComparison function) 1x1
    %         string
    %       - met_idx: array 1xn, every element in the array is a double
    %         1xn with the positions of the mets to be used in the fluxsum
    %         comptutation, the positions need to be the position of that
    %         metabolite in the choosen reference model
    %       - rxn_idx: same as met_idx
    %       - rxn_set_labels: names of the defined sets, can be choosen
    %         freely
    %       - plot_type: to visualize a heatmap (fluxsum overall given in
    %         one rnx_id set) or violin (fluxsum for all metabolites taking
    %         part in the specified rxns)
    %       - exclude_coenzymes: option to exclude all the metabolite that
    %         have a connectivity of over 135 within the reference network
    %         (those are most likely coenzymes -> like h2o)
    %       - ignore_compartment: allows to ignore the compartment and
    %         build an overall fluxsum per metabolite (only relevant for
    %         violin plots)
    % Output: 
    %       - fluxsum_sets: fluxsum overall metabolites per sample as a
    %       array with matrices stored within
    arguments
        project struct
        comparison_name (1,1) string
        met_idx (1,:) cell {mustBeColumnVector} =[]
        rxn_idx (1,:) cell {mustBeColumnVector} =[]
        rxn_set_labels (1,:) string = []
        plot_type  {mustBeMember(plot_type, ["violin","heatmap_model", "heatmap_sample", "heatmap_sample_all_mets"])} =["violin"] 
        exclude_coenzymes (1,1) logical = true
        ignore_compartment (1,1) logical = true
    end
    
    % sets the order of the rxns in order to be able to bring the samples
    % from all the models in the same order
    reference = project.comparisons.(comparison_name).reference_model;

    % by default if no input is given the fluxsum for all metabolites overa
    % ll reaction is computed
    if isempty(rxn_idx)
        rxn_idx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    if isempty(met_idx)
        met_idx = find(ones(length(project.models.(reference).model.mets),1));
    end
    
    % in order to prevent coenzymes from inflating the fluxsum of specific
    % pathways all metabolites which have a connectivity higher than 135 in
    % the choosen reference model will be excluded from the computation of
    % the fluxsum and not be visualized in the violinplot
    if exclude_coenzymes
        S = project.models.(reference).model.S;
        rxn_count_per_met = sum(S ~= 0, 2);  % sum over columns = number of reactions each metabolite participates in
        idx_coenzymes = find(rxn_count_per_met >135);

        % histogram(rxn_count, 1000)
        % xlim([0,50])
        % xlabel('Number of reactions per metabolite')
        % ylabel('Count of metabolites')
        
        % define coenzymes by connectivity in the network, highly connected
        % metabolites are probably coenzymes 
        % idx_coenzymes = find(rxn_count_per_met > quantile(rxn_count_per_met,0.994));

        met_idx = setdiff(met_idx,idx_coenzymes);
    end
    
    % depending on the choosen plot type different functions are executed
    if plot_type == "violin"
            fluxsum_sets = get_violin_plots(project,comparison_name,met_idx, rxn_idx,rxn_set_labels,ignore_compartment);
    else
            fluxsum_sets = get_comparison_heatmap(project,comparison_name,met_idx, rxn_idx,rxn_set_labels,plot_type);
    end

    % when met_idx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites
end

function fluxsum_sets = get_comparison_heatmap(project,comparison_name,met_idx,rxn_idx,rxn_set_labels,type)
    % This function visualizes the mean fluxsum over the specified rxn_id
    % sets. By getting the fluxsum for metabolites that are participating
    % in the specified sets + are part of the met_idx and then computing
    % the mean over each set. Using this function a comparison over models
    % and different defined subsystems can be visualized. 
    % INPUT:
    %       - project: project structure after running modelsComparison
    %         function 
    %       - comparison_name: name of the comparison you want to visualize
    %         (the second output of the modelsComparison function) 1x1
    %         string
    %       - met_idx: with the positions of the mets to be used in the fluxsum
    %         comptutation, the positions need to be the position of that
    %         metabolite in the choosen reference model
    %       - rxn_idx: sets of rxns to be visualized in the heatmap
    %       - rxn_set_labels: names of the defined sets, can be choosen
    %         freely
    %       - type: specifies whether a fluxsum per rxn set should be
    %         visualized per samples or the average over all samples from one
    %         model: values : "heatmap_sample", "heatmap_sample_all_mets" or "heatmap_model" default is "heatmap_model"
    % Output: 
    %       - fluxsum_sets: fluxsum overall metabolites per sample as a
    %         array with matrices stored within
    arguments
        project struct
        comparison_name (1,1) string
        met_idx  
        rxn_idx
        rxn_set_labels (1,:) string 
        type {mustBeMember(type, ["heatmap_model", "heatmap_sample","heatmap_sample_all_mets"])} =["heatmap_model"] 
    end
    reference = project.comparisons.(comparison_name).reference_model;
    if isempty(rxn_idx)
        rxn_idx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    
    if isempty(met_idx)
        met_idx = find(ones(length(project.models.(reference).model.mets),1));
    end
    

    fluxsum_sets = get_fluxsum(project,comparison_name,met_idx,rxn_idx);
    samples_cat = cellstr(project.comparisons.(comparison_name).sample_model_labels);


    heatmap_data = zeros(length(fluxsum_sets), length(unique(samples_cat)));
    heatmap_data_all_samples = zeros(length(fluxsum_sets), length(samples_cat));

    % when met_idx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites
    heatmap_data_all_samples_all_mets = {};
    for subsystem = 1:numel(fluxsum_sets)
        data = fluxsum_sets{subsystem};
        title_fig = rxn_set_labels(subsystem);

        met_names = project.models.(reference).model.mets(met_idx);
    
        met_names = met_names(find(any(data ~= 0,2)));
        data = data(find(any(data ~= 0,2)),:);
        
    
        %
        samples_cat = categorical(samples_cat);  % now KO/PLV/WT are categories
        groups = categories(samples_cat);        % {"KO","PLV","WT"}
        nGroups = numel(groups);
        [nMet, nSamples] = size(data);
        data_grouped = cell(1,nGroups);
    
        for g = 1:nGroups
            idx = samples_cat == groups(g);   % logical index for this group
            data_grouped{g} = data(:, idx);   % all metabolites, only this group
        end
        heatmap_data(subsystem,:) = cellfun(@(x) mean(x(:)), data_grouped);
        heatmap_data_all_samples(subsystem,:) = cell2mat(cellfun(@(x) mean(x,1), data_grouped, 'UniformOutput', false));
        heatmap_data_all_samples_all_mets{end +1} = cell2mat(data_grouped);
    end

    if type == "heatmap_model"
        % z-scaling for the heatmap in order to make the differences between
        % samples for one pathway more visible!
        scaled_data = zscore(heatmap_data')';
    
        imagesc(scaled_data)
        
        cmap = get_color_pallette();
        h = colorbar;  
        caxis([-max([max(scaled_data(:)),abs(min(scaled_data(:)))]) max([max(scaled_data(:)),abs(min(scaled_data(:)))])])   % colors scaled from -2 (min) to 2 (max)

        ylabel(h, 'Scaled average fluxsum average value over rxn set', 'FontSize', 18)        % Set title/label of colorbar
        axis equal tight            % Make cells square and remove extra space
    
        title("Fluxsum per rxn set overall metabolites")    % grayscale
        % Set x-axis and y-axis labels
        set(gca, 'XTick', 1:length(unique(samples_cat)), 'XTickLabel', unique(samples_cat), ...
             'YTick', 1:length(rxn_set_labels), 'YTickLabel', rxn_set_labels)
        xtickangle(45)
        ax = gca;
        ax.FontSize = 18;  
        xlabel('Model', 'FontSize', 18)       
        ylabel('Reaction set', 'FontSize', 18)    
       
        [nRows, nCols] = size(heatmap_data);
        
    
        % Loop over every cell and place the absolute number from pathway_counts
        for i = 1:nRows
            for j = 1:nCols
                % You want absolute numbers, not relative counts, so use pathway_counts
                % (or multiply relative_counts by reference_model if needed)
                value =  heatmap_data(i,j); % +1 because first column is reference_model
                % Place text at the center of the tile
                text(j, i, num2str(value), ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'Color','k', ...          % black text
                    'FontSize',18)
            end
        end
        
        hold off
    elseif type == "heatmap_sample"

        figure
        scaled_data = zscore(heatmap_data_all_samples')';
    
        imagesc(scaled_data)

        cmap = get_color_pallette();
       
        h = colorbar;  
        caxis([-max([max(scaled_data(:)),abs(min(scaled_data(:)))]) max([max(scaled_data(:)),abs(min(scaled_data(:)))])])   % colors scaled from -2 (min) to 2 (max)

        ylabel(h, 'Scaled average fluxsum value per sample', 'FontSize', 18)        % Set title/label of colorbar
        title("Fluxsum per rxn set overall metabolites")    % grayscale
        % Set x-axis and y-axis labels
        [sample_count,~] = hist(samples_cat);
        xtickposition = ((1:length(unique(samples_cat))) .* (sample_count)) - sample_count/2;
    
        set(gca, 'XTick', xtickposition, 'XTickLabel', unique(samples_cat), ...
             'YTick', 1:length(rxn_set_labels), 'YTickLabel', rxn_set_labels)
        xtickangle(45)
        ax = gca;
        ax.FontSize = 18;  
        xlabel('Model', 'FontSize', 18)       
        ylabel('Reaction set', 'FontSize', 18)
    else
    

    % z-scaling for the heatmap in order to make the differences between
        % samples for one pathway more visible!
        scaled_data = zscore(cell2mat(heatmap_data_all_samples_all_mets')')';
        
    
        imagesc(scaled_data)
        
        cmap = get_color_pallette();
        h = colorbar;  
        min_axis = quantile(scaled_data(:),0.001);
        max_axis = quantile(scaled_data(:),0.999);
        axis_limit = max([abs(min_axis), abs(max_axis)]);
        caxis([-axis_limit,axis_limit])   

        ylabel(h, 'Scaled average fluxsum average value in rxn set', 'FontSize', 18)        % Set title/label of colorbar
        
        title("Fluxsum per rxn set overall metabolites")    % grayscale
        % Set x-axis and y-axis labels
        [sample_count,~] = hist(samples_cat);
        xtickposition = ((1:length(unique(samples_cat))) .* (sample_count)) - sample_count/2;
    
        count_mets_per_set = cell2mat(arrayfun(@(x)size(x{:},1),heatmap_data_all_samples_all_mets,'UniformOutput',false))'
        ytickposition = cumsum(count_mets_per_set) - round(count_mets_per_set/2);
    
        set(gca, 'XTick', xtickposition, 'XTickLabel', unique(samples_cat), ...
             'YTick', ytickposition, 'YTickLabel', rxn_set_labels)
        xtickangle(45)
        ax = gca;
        ax.FontSize = 18;  
        xlabel('Model', 'FontSize', 18)       
        ylabel('Reaction set', 'FontSize', 18)    
    end 

end

function fluxsum_sets = get_violin_plots(project,comparison_name,met_idx,rxn_idx,rxn_set_labels, ignore_compartment,model_name_flux_boundaries)
    % This function visualizes the fluxsum distribution per metabolite and model 
    % over the specified rxn_id sets into a violin plot.
    % By getting the fluxsum for metabolites that are participating
    % in the specified sets + are part of the met_idx.
    % Using this function a comparison over models for the fluxsum of
    % different metabolites can be visualized.  
    % INPUT:
    %       - project: project structure after running modelsComparison
    %         function 
    %       - comparison_name: name of the comparison you want to visualize
    %         (the second output of the modelsComparison function) 1x1
    %         string
    %       - met_idx: with the positions of the mets to be used in the fluxsum
    %         comptutation, the positions need to be the position of that
    %         metabolite in the choosen reference model
    %       - rxn_idx: rxn sets, for each set a separate figure is created
    %         with all the metabolite fluxsum values for all the metabolite
    %         participating in the specified rxn_idx
    %       - rxn_set_labels: names of the defined sets, can be choosen
    %         freely
    %       - ignore_compartment: allows to compute the fluxsum for a
    %         metabolite ignoring the cellular location
    %       - model_name_flux_boundaries: model that saves the
    %         constraints set when constructing the models
    % Output: 
    %       - fluxsum_sets: fluxsum overall metabolites per sample as a
    %         array with matrices stored within
    arguments
        project struct
        comparison_name (1,1) string
        met_idx
        rxn_idx
        rxn_set_labels (1,:)
        ignore_compartment (1,1) logical = true
        model_name_flux_boundaries (1,1) string = "consistent_medium_constrained_model"
    end

    reference = project.comparisons.(comparison_name).reference_model;
    reference_model = project.models.(reference).model;
    modelNames = project.comparisons.(comparison_name).modelNames;
    models = rmfield(project.models,setdiff(fieldnames(project.models),...
                                            modelNames));
    if isempty(rxn_idx)
        rxn_idx = find(ones(reference_model.rxns),1);
    end
    
    if isempty(met_idx)
        met_idx = find(ones(length(reference_model.mets),1));
    end
    

    fluxsum_sets = get_fluxsum(project,comparison_name,met_idx,rxn_idx);

    
    % for every specified set in rxn_idx create one figure with all the
    % fluxsums of the metabolites participating in the rxns + being part of
    % met_idx
    for subsystem = 1:numel(fluxsum_sets)
        data = fluxsum_sets{subsystem};
        title_fig = rxn_set_labels(subsystem);
        
        met_names = reference_model.mets(met_idx); % filter for mets in met_idx 
        
        zero_fluxsum_metabolites = find(any(data ~= 0,2)); % filter for metabolites which have zero fluxsum overall samples, overall models
        data = data(zero_fluxsum_metabolites,:);
        met_names = met_names(zero_fluxsum_metabolites);
        samples_cat = cellstr(project.comparisons.(comparison_name).sample_model_labels);
        % remove the compartment specification to compute the totall
        % fluxsum of a metabolite, and add them up overall compartments
        if ignore_compartment
            met_names_without_compartment = regexprep(met_names,"\[.\]$", "");

            % Get unique specifications and grouping indices
            [unique_mets, ~, idx] = unique(met_names_without_compartment);
            
            % Prepare output
            numSpecs  = numel(unique_mets);
            numCols   = size(data, 2); 
            aggMatrix = zeros(numSpecs, numCols);
            
            % Sum rows by group
            for i = 1:numSpecs
                agg_data(i, :) = sum(data(idx == i, :), 1);
            end

            data = agg_data;
            met_names = unique_mets;

        end
        
        
        samples_cat = categorical(samples_cat);  
        groups = categories(samples_cat);       
        nGroups = numel(groups);
        [nMet, nSamples] = size(data);
        data_grouped = cell(1,nGroups);
        
        % put the samples from different models in different matrices
        for g = 1:nGroups
            idx = samples_cat == groups(g);   
            data_grouped{g} = data(:, idx);   
        end
        
        % put each metabolite in a different vector
        for m = 1:nMet
            for g = 1:nGroups
                data_for_plot{m,g} = data_grouped{g}(m,:);  % 1 × nSamples_in_group
            end
        end
        % now we have a separate vector for every metabolite for every
        % model
        % now we loop over all metabolites to visualize on violin per model
        
        % Compute rows and columns to make layout roughly square
        nRows = ceil(sqrt(nMet));
        nCols = ceil(nMet / nRows);
        
        % Create figure with quadratic tiled layout
        figure
        tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact')
        
        for m = 1:nMet
            ax = nexttile;    % get the axes handle
            hold(ax, 'on')
        
            dat = data_for_plot(m,:);        
            dat = vertcat(dat{:})';         
            
            columns_to_keep = find(~all(dat == 0));
            violinplot(dat(:,columns_to_keep), groups(columns_to_keep),'ShowData', false);    % leave out 'ShowData' for compatibility
       
            title(met_names{m}, 'Interpreter','none');
            ylabel(ax, 'Flux Value', 'FontSize', 18);
            title(ax, met_names{m}, 'Interpreter','none', 'FontSize', 18);
            ax.FontSize = 18;     
        end
    
        sgtitle("Metabolite Fluxsum in subSystem: " + title_fig,'FontSize',20)

        % get Table with rxn formula, discretization status, concentration
        % in the medium etc
        T_display = get_rxn_overview_table(project, comparison_name,models,reference_model, met_names,model_name_flux_boundaries,modelNames,reference);
        
        T_display = addvars(T_display, string(T_display.Properties.RowNames), ...
                            'Before', 1, ...
                            'NewVariableNames', "Reaction");
        
        % Create UI figure
        fig = uifigure('Name', "Metabolite Fluxsum in subSystem: " + title_fig, ...
                       'Position',[100 100 1200 600]);
        
        % Create table filling the entire figure
        tbl = uitable(fig, ...
            'Data', T_display{:,:}, ...
            'ColumnName', T_display.Properties.VariableNames, ...
            'Units','normalized', ...
            'Position',[0 0 1 1], ...
            'FontSize',18, ...
            'ColumnWidth','auto');

        

    end

end

function T_display = get_rxn_overview_table(project, comparison_name,models,reference_model, met_names,model_name_flux_boundaries,modelNames,reference )
        
        ref = fieldnames(models);
        ref = models.(ref{1,1}).settings.medium;
        
        % get rxn_idx
        [rxn_names] = findRxnsFromMets(reference_model,string(met_names));
        rxn_ids = find(matches(reference_model.rxns,rxn_names));
         
        % get rxn_idx formulas + filter out rxns that have zero in all
        % samples
        samples = project.comparisons.(comparison_name).ordered_samples;
        zero_rxns = find(sum(samples == 0,2) == size(samples,2));
        rxn_ids = setdiff(rxn_ids,zero_rxns);
        rxn_names = reference_model.rxns(rxn_ids);
        
        % get lower and upper bound for every rxns 
        ordered_lb = getOrderedFeatureMatrix(project,model_name_flux_boundaries,...
                                             "rxns",reference,"model.lb");
        ordered_ub = getOrderedFeatureMatrix(project,model_name_flux_boundaries,...
                                             "rxns",reference,"model.ub");
        ordered_mapping_rxn_matrix = getOrderedFeatureMatrix(project,modelNames,...
                                                             "rxns",reference,"mappedDiscRxns");
    
    
        % check with samplings are not zero 

        medium_constrained = ismember(rxn_names, ref.medium_composition.ExRxns_Recon3D) | ...
                         ismember(rxn_names, ref.manual_set_boundaries.unwanted_export) | ...
                         ismember(rxn_names, ref.manual_set_boundaries.unwanted_import);

        ordered_ub = ordered_ub(rxn_ids,:);
        ordered_lb = ordered_lb(rxn_ids,:);
        ordered_mapping_rxn_matrix = ordered_mapping_rxn_matrix(rxn_ids,:);
        rxn_formulas = string(printRxnFormula(reference_model));
        rxn_formulas = rxn_formulas(rxn_ids);
  
        % get rxn gene rules to add to the table
        symbol_gpr_rules = string(cellfun(@(rxnName)get_rxn_symbol_rule(project.models.(reference),...
                                                   rxnName),string(rxn_names),'UniformOutput', false));

        T = table(rxn_formulas, medium_constrained,...
            join(string(ordered_mapping_rxn_matrix), "|", 2),symbol_gpr_rules,...
                  'VariableNames', ["Reaction Formula","medium constrained",...
                                    join(string(project.comparisons.(comparison_name).modelNames),"_"),...
                                    "symbol gpr rules"], ...
                  'RowNames',rxn_names);
        T = T(flip(string(T.Properties.RowNames)),:);
        T.lb = flip(ordered_lb);
        T.ub = flip(ordered_ub);

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



