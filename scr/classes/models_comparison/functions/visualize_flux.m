function flux_sets = visualize_flux(project,comparison_name,met_idx,rxn_idx,rxn_set_labels,threshold)
    arguments
        project
        comparison_name
        met_idx =[]
        rxn_idx =[]
        rxn_set_labels = []
        threshold {mustBeMember(threshold, ["all","positive","negative"])} =["all"] 
    end
    %     plot_type  {mustBeMember(plot_type, ["violin","heatmap"])} =["violin"] 
    %     exclude_coenzymes (1,1) logical = true
    %     ignore_compartment (1,1) logical = true
    % end
    reference = project.comparisons.(comparison_name).reference_model;
    if isempty(rxn_idx)
        rxn_idx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    
    if isempty(met_idx)
        met_idx = find(ones(length(project.models.(reference).model.mets),1));
    end

    flux_sets = get_violin_plots_flux(project,comparison_name,met_idx, rxn_idx,rxn_set_labels,threshold);

end

function flux_sets = get_violin_plots_flux(project,comparison_name,met_idx,rxn_idx,rxn_set_labels, threshold)
    arguments
        project
        comparison_name
        met_idx =[]
        rxn_idx =[]
        rxn_set_labels = []
        threshold {mustBeMember(threshold, ["all","positive","negative"])} =["positive"] 
    end

    reference = project.comparisons.(comparison_name).reference_model;
    reference_model = project.models.(reference).model;
    modelNames = project.comparisons.(comparison_name).modelNames;
    models = rmfield(project.models,setdiff(fieldnames(project.models),...
                                            modelNames));
    if isempty(rxn_idx)
        rxn_idx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    
    if isempty(met_idx)
        met_idx = find(ones(length(project.models.(reference).model.mets),1));
    end
    

    flux_sets = get_flux(project,comparison_name,met_idx,rxn_idx);
    
    % when met_idx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites

    

    for subsystem = 1:numel(flux_sets)
        data = flux_sets{subsystem};
        title_fig = rxn_set_labels(subsystem);

        rxns_names = project.models.(reference).model.rxns(rxn_idx{subsystem});
        

        zero_fluxsum_rxns = find(any(data ~= 0,2));
        data = data(zero_fluxsum_rxns,:);
        rxns_names = rxns_names(zero_fluxsum_rxns);

        if threshold ~= "all"
            Q = zeros(128, 4);
            for i = 1:size(data,1)
                Q(i,:) = quantile(data(i,:), [0.20, 0.40,0.60,0.80]);
            end
            % compute quantiles -> check that at least 3  out of the 4
            % quantiles are either positive or negative

            if threshold == "positive"
                idx = find(sum(Q > 0,2) >2);
            elseif threshold == "negative"
                idx = find(sum(Q < 0,2) >2);
            end

            data = data(idx,:);
            rxns_names = rxns_names(idx);
        end

        
        samples_cat = cellstr(project.comparisons.(comparison_name).sample_model_labels);
        
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
        nCols = ceil(nMet / nRows);
        
        %% Create figure with quadratic tiled layout
        figure
        tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact')
        
        for m = 1:nMet
            ax = nexttile;
            hold on
        
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
            ylabel(ax, 'Flux Value', 'FontSize', 18)
            title(ax, rxns_names{m}, 'Interpreter','none', 'FontSize', 18)
            ax.FontSize = 18;     % sets tick labels
        end
    
        sgtitle("Rxns Flux in subSystem: " + title_fig,'FontSize',20)


        %% table with rxnforumlas etc 

        ref = fieldnames(models);
        ref = models.(ref{1,1}).settings.medium;
        %media_for_models = structfun(@(x) x.settings.medium, models);
        %medium_is_equal_between_models = all(arrayfun(@(x) isequaln(x, ref), media_for_models));

        rxn_ids = find(matches(string(reference_model.rxns),rxns_names));
        

        samples = project.comparisons.(comparison_name).ordered_samples;
        zero_rxns = find(sum(samples == 0,2) == size(samples,2));
        rxn_ids = setdiff(rxn_ids,zero_rxns);
        rxns_names = reference_model.rxns(rxn_ids);
        

        ordered_lb = getOrderedFeatureMatrix(project,"consistent_medium_constrained_model",...
                                             "rxns",reference,"model.lb");
        ordered_ub = getOrderedFeatureMatrix(project,"consistent_medium_constrained_model",...
                                             "rxns",reference,"model.ub");
        ordered_mapping_rxn_matrix = getOrderedFeatureMatrix(project,modelNames,...
                                                             "rxns",reference,"mappedDiscRxns");
    
    
        % check with samplings are not zero 

        medium_constrained = ismember(rxns_names, ref.medium_composition.ExRxns_Recon3D) | ...
                         ismember(rxns_names, ref.manual_set_boundaries.unwanted_export) | ...
                         ismember(rxns_names, ref.manual_set_boundaries.unwanted_import);

        ordered_ub = ordered_ub(rxn_ids,:);
        ordered_lb = ordered_lb(rxn_ids,:);
        ordered_mapping_rxn_matrix = ordered_mapping_rxn_matrix(rxn_ids,:);
        rxn_formulas = string(printRxnFormula(reference_model));
        rxn_formulas = rxn_formulas(rxn_ids);
  
        % get rxn gene rules to add to the table
        symbol_gpr_rules = string(cellfun(@(rxnName)get_rxn_symbol_rule(project.models.(reference),...
                                                   rxnName),string(rxns_names),'UniformOutput', false));

        T = table(rxn_formulas, medium_constrained,...
            join(string(ordered_mapping_rxn_matrix), "|", 2),symbol_gpr_rules,...
                  'VariableNames', ["Reaction Formula","medium constrained",...
                                    join(string(project.comparisons.(comparison_name).modelNames),"_"),...
                                    "symbol gpr rules"], ...
                  'RowNames',rxns_names);
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

function flux_cell = get_flux(project,comparison_name,met_idx,rxn_idx)
    arguments
        project
        comparison_name
        met_idx =[]
        rxn_idx =[]
    end
    reference = project.comparisons.(comparison_name).reference_model;
    if isempty(rxn_idx)
        rxn_idx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    
    if isempty(met_idx)
        met_idx = find(ones(length(project.models.(reference).model.mets),1));
    end
    

    samples = project.comparisons.(comparison_name).ordered_samples;
    
    flux_cell = cellfun(@(x) samples(x,:), ...
                           rxn_idx,'UniformOutput', false);

end

