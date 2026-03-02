function fluxsum_sets = visualize_fluxsum(project,comparison_name,met_idx,rxn_idx,rxn_set_labels,plot_type,exclude_coenzymes,ignore_compartment)
    arguments
        project
        comparison_name
        met_idx =[]
        rxn_idx =[]
        rxn_set_labels = []
        plot_type  {mustBeMember(plot_type, ["violin","heatmap"])} =["violin"] 
        exclude_coenzymes (1,1) logical = true
        ignore_compartment (1,1) logical = true
    end
    reference = project.comparisons.(comparison_name).reference_model;

    if isempty(rxn_idx)
        rxn_idx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    
    if isempty(met_idx)
        met_idx = find(ones(length(project.models.(reference).model.mets),1));
    end

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
    

    if plot_type == "violin"
            fluxsum_sets = get_violin_plots(project,comparison_name,met_idx, rxn_idx,rxn_set_labels,ignore_compartment);
    else
            fluxsum_sets = get_comparison_heatmap(project,comparison_name,met_idx, rxn_idx,rxn_set_labels);
    end

    
    
    % when met_idx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites
end

function fluxsum_sets = get_comparison_heatmap(project,comparison_name,met_idx,rxn_idx,rxn_set_labels)
    arguments
        project
        comparison_name
        met_idx =[]
        rxn_idx =[]
        rxn_set_labels = []
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
    
    % when met_idx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites
    for subsystem = 1:numel(fluxsum_sets)
        data = fluxsum_sets{subsystem};
        title_fig = rxn_set_labels(subsystem);

        met_names = project.models.(reference).model.mets(met_idx);
    
        data = data(find(any(data ~= 0,2)),:);
        met_names = met_names(find(any(data ~= 0,2)));
    
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
        
        
    end

    % Optional row/column labels
    row_labels = rxn_set_labels;
    col_labels = groups;
    
    % Create heatmap
    figure;
    hold off;
    h = heatmap(col_labels, regexprep(row_labels,"_"," "), heatmap_data);
    
    % Customize
    h.Title = 'Flux Sum Heatmap';
    h.XLabel = 'Model';
    h.YLabel = 'pathway';
    h.Colormap = parula;        % or 'hot', 'jet', etc.
    h.ColorLimits = [min(data(:)) max(data(:))];  % optional
    

end

function fluxsum_sets = get_violin_plots(project,comparison_name,met_idx,rxn_idx,rxn_set_labels, ignore_compartment)
    arguments
        project
        comparison_name
        met_idx =[]
        rxn_idx =[]
        rxn_set_labels = []
        ignore_compartment (1,1) logical = true
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
    
    % when met_idx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites

    

    for subsystem = 1:numel(fluxsum_sets)
        data = fluxsum_sets{subsystem};
        title_fig = rxn_set_labels(subsystem);

        met_names = reference_model.mets(met_idx);
        
        zero_fluxsum_metabolites = find(any(data ~= 0,2));
        data = data(zero_fluxsum_metabolites,:);
        met_names = met_names(zero_fluxsum_metabolites);
        samples_cat = cellstr(project.comparisons.(comparison_name).sample_model_labels);
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
    
        
        for m = 1:nMet
            for g = 1:nGroups
                data_for_plot{m,g} = data_grouped{g}(m,:);  % 1 × nSamples_in_group
            end
        end
        
        % Compute rows and columns to make layout roughly square
        nRows = ceil(sqrt(nMet));
        nCols = ceil(nMet / nRows);
        
        %% Create figure with quadratic tiled layout
        figure
        tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact')
        
        for m = 1:nMet
            ax = nexttile;    % get the axes handle
            hold(ax, 'on')
        
            % Step 1: prepare data for violinplot
            dat = data_for_plot(m,:);        % 1 x nGroups cell
            dat = vertcat(dat{:})';          % column vector of all samples
        
            % Step 2: create group vector for violinplot
            group_vector = repelem(groups, cellfun(@numel, data_for_plot(m,:)))';
        
            % Step 3: plot
            columns_to_keep = find(~all(dat == 0));
            obj = violinplot(dat(:,columns_to_keep), groups(columns_to_keep),'ShowData', false);    % leave out 'ShowData' for compatibility
        
            % Step 4: labels
            % xticks(1:numel(groups))
            % xticklabels(groups)
            
            title(met_names{m}, 'Interpreter','none')
            ylabel(ax, 'Flux Value', 'FontSize', 18)
            title(ax, met_names{m}, 'Interpreter','none', 'FontSize', 18)
            ax.FontSize = 18;     % sets tick labels
        end
    
        sgtitle("Metabolite Fluxsum in subSystem: " + title_fig,'FontSize',20)

        %% table with rxnforumlas etc 

        ref = fieldnames(models);
        ref = models.(ref{1,1}).settings.medium;
        %media_for_models = structfun(@(x) x.settings.medium, models);
        %medium_is_equal_between_models = all(arrayfun(@(x) isequaln(x, ref), media_for_models));
    

        
        
        [rxn_names] = findRxnsFromMets(reference_model,string(met_names));
        rxn_ids = find(matches(reference_model.rxns,rxn_names));
        

        samples = project.comparisons.(comparison_name).ordered_samples;
        zero_rxns = find(sum(samples == 0,2) == size(samples,2));
        rxn_ids = setdiff(rxn_ids,zero_rxns);
        rxn_names = reference_model.rxns(rxn_ids);
        

        ordered_lb = getOrderedFeatureMatrix(project,"consistent_medium_constrained_model",...
                                             "rxns",reference,"model.lb");
        ordered_ub = getOrderedFeatureMatrix(project,"consistent_medium_constrained_model",...
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





% function project = compute_flux_sum(project,list_model_names,reference_model,compute_based_on_incoming_flux)
%             arguments
%                 project
%                 list_model_names
%                 reference_model
%                 compute_based_on_incoming_flux =0
%             end
% 
% 
%     for mod = list_model_names
%         model = project.models.(mod);
%         project.models.(mod).analysis.sampling.fluxsum = get_model_fluxsum(model,model.analysis.sampling.samples);
%     end
% 
% 
% end