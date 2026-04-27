function [fluxSets,figs] = visualizeFlux(project,comparisonName,metIdx,rxnIdx,rxnSetLabels,threshold,plotVisible)
    arguments
        project
        comparisonName
        metIdx =[]
        rxnIdx =[]
        rxnSetLabels = []
        threshold {mustBeMember(threshold, ["all","positive","negative"])} =["all"] 
        plotVisible ="on"
    end
    %     plot_type  {mustBeMember(plot_type, ["violin","heatmap"])} =["violin"] 
    %     exclude_coenzymes (1,1) logical = true
    %     ignore_compartment (1,1) logical = true
    % end
    reference = project.comparisons.(comparisonName).referenceModel;
    if isempty(rxnIdx)
        rxnIdx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    
    if isempty(metIdx)
        metIdx = find(ones(length(project.models.(reference).model.mets),1));
    end

    [fluxSets,figs] = getViolinPlotsFlux(project,comparisonName,metIdx, rxnIdx,rxnSetLabels,threshold,plotVisible);

end

function [fluxSets,figs] = getViolinPlotsFlux(project,comparisonName,metIdx,rxnsIdx,rxnSetLabels, threshold,plotVisible)
    arguments
        project
        comparisonName
        metIdx =[]
        rxnsIdx =[]
        rxnSetLabels = []
        threshold {mustBeMember(threshold, ["all","positive","negative"])} =["positive"] 
        plotVisible ="on"
    end

    reference = project.comparisons.(comparisonName).referenceModel;
    referenceModel = project.models.(reference).model;
    modelNames = project.comparisons.(comparisonName).modelNames;
    models = rmfield(project.models,setdiff(fieldnames(project.models),...
                                            modelNames));
    if isempty(rxnsIdx)
        rxnsIdx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    
    if isempty(metIdx)
        metIdx = find(ones(length(project.models.(reference).model.mets),1));
    end
    

    fluxSets = getFlux(project,comparisonName,metIdx,rxnsIdx);
    
    % when met_idx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites

    figs = struct();

    for subsystem = 1:numel(fluxSets)
        data = fluxSets{subsystem};
        title_fig = rxnSetLabels(subsystem);
        plot_name = replace(title_fig, ["_", "-", "/"], "");

        rxns_names = project.models.(reference).model.rxns(rxnsIdx{subsystem});
        

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

        
        samples_cat = project.comparisons.(comparisonName).sampleModelLabels;
        groups = unique(samples_cat, 'stable');       
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
        
        maxPlotsPerFig = 12;
        nFigs = ceil(nMet / maxPlotsPerFig);
        
        figStruct = struct();

        for f = 1:nFigs
            
            fig = figure('Color','w','Position',[100 100 800 800],...
                                       'Visible',plotVisible);
            
            t = tiledlayout(3,4);
            title(t,"Flux for reactions in : " + title_fig, ...
                  'FontSize', 20, 'FontWeight','bold', 'Interpreter','none');
            
            % Store using dynamic field name
            fieldName = sprintf('fig%d', f);
            figStruct.(fieldName) = fig;
            
            startIdx = (f-1)*maxPlotsPerFig + 1;
            endIdx   = min(f*maxPlotsPerFig, nMet);
            
            for m = startIdx:endIdx
                
                ax = nexttile(t);
                hold(ax, 'on')
                
                dat = data_for_plot(m,:);
                % in order to be able to put it all in one matrix we need
                % to make the sample count the same length, so we add NaNs
                
                maxLen = max(cellfun(@numel, dat));
                dat = cell2mat(cellfun(@(x) [x, nan(1, maxLen-numel(x))], dat, 'UniformOutput', false));
                dat = reshape(dat, maxLen, []);
                
                
                columns_to_keep = find(~all(dat == 0 | isnan(dat)));
                evalc('violinplot(dat(:,columns_to_keep), groups(columns_to_keep), "ShowData", false);');
                
                ylabel(ax, 'Flux Value', 'FontSize', 18);
                title(ax, rxns_names{m}, 'Interpreter','none', 'FontSize', 14);
                
                ax.FontSize = 18;
            end
        end
        
        if plotVisible == "on"
    

        %table with rxnforumlas etc 

        ref = fieldnames(models);
        ref = models.(ref{1,1}).settings.medium;
        %media_for_models = structfun(@(x) x.settings.medium, models);
        %medium_is_equal_between_models = all(arrayfun(@(x) isequaln(x, ref), media_for_models));

        rxn_ids = find(matches(string(referenceModel.rxns),rxns_names));
        

        samples = project.comparisons.(comparisonName).orderedSamples;
        zero_rxns = find(sum(samples == 0,2) == size(samples,2));
        rxn_ids = setdiff(rxn_ids,zero_rxns);
        rxns_names = referenceModel.rxns(rxn_ids);
        

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
       
        rxn_abbr = referenceModel.rxns(rxn_ids);
        rxn_formulas = string(printRxnFormula(referenceModel,rxn_abbr,false));
  
        % get rxn gene rules to add to the tabler
        symbol_gpr_rules = string(cellfun(@(rxnName)getRxnSymbolRule(project.models.(reference),...
                                                   rxnName),string(rxns_names),'UniformOutput', false));

        T = table(rxn_formulas, medium_constrained,...
            join(string(ordered_mapping_rxn_matrix), "|", 2),symbol_gpr_rules,...
                  'VariableNames', ["Reaction Formula","medium constrained",...
                                    join(string(project.comparisons.(comparisonName).modelNames),"_"),...
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
        plot_name = replace(plot_name, ["_", "-", "/", " "], "");
        figs.("violinFlux" + plot_name) = figStruct;

    end

end

function fluxCell = getFlux(project,comparisonName,metIdx,rxnIdx)
    arguments
        project
        comparisonName
        metIdx =[]
        rxnIdx =[]
    end
    reference = project.comparisons.(comparisonName).referenceModel;
    if isempty(rxnIdx)
        rxnIdx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    
    if isempty(metIdx)
        metIdx = find(ones(length(project.models.(reference).model.mets),1));
    end
    

    samples = project.comparisons.(comparisonName).orderedSamples;
    
    fluxCell = cellfun(@(x) samples(x,:), ...
                           rxnIdx,'UniformOutput', false);

end

