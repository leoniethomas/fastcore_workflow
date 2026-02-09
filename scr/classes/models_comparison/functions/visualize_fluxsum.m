function visualize_fluxsum(project,comparison_name,met_idx,rxn_idx,rxn_set_labels)
    arguments
        project
        comparison_name
        met_idx =[]
        rxn_idx =[]
        rxn_set_labels = [];
    end
    reference = project.comparisons.(comparison_name).reference_model;
    if isempty(rxn_idx)
        rxn_idx = find(ones(length(project.models.(reference).model.rxns),1));
    end
    
    if isempty(met_idx)
        met_idx = find(ones(length(project.models.(reference).model.mets),1));
    end
    

    fluxsum_sets = get_fluxsum(project,comparison_name,met_idx,rxn_idx);
    
    % when met_idx is empyt, or over a specific number of mets -> over 50
    % then only display the top metabolites

    

    for subsystem = 1:numel(fluxsum_sets)
        data = fluxsum_sets{subsystem};
        title_fig = rxn_set_labels(subsystem);

        met_names = project.models.(reference).model.mets(met_idx);
    
        data = data(find(any(data ~= 0,2)),:);
        met_names = met_names(find(any(data ~= 0,2)));
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
        
        %% Create figure with quadratic tiled layout
        figure
        tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact')
        
        for m = 1:nMet
            nexttile
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
            ylabel('Value')
            title(met_names{m}, 'Interpreter','none')
        end
    
        sgtitle("Metabolite Fluxsum in subSystem: " + title_fig,'FontSize',20)
    end

end

function fluxsum_cell = get_fluxsum(project,comparison_name,met_idx,rxn_idx)
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
    full_S = project.models.(reference).model.S;
    
    fluxsum_cell = cellfun(@(x) get_model_fluxsum(full_S(:,x), samples(x,:)), ...
                           rxn_idx,'UniformOutput', false);

    fluxsum_cell = cellfun(@(f) f(met_idx,:), ...
                           fluxsum_cell, 'UniformOutput', false);
    

end


function fluxsum = get_model_fluxsum(S_matrix,samples,flux_summed_up)
    arguments
        S_matrix 
        samples 
        flux_summed_up {mustBeMember(flux_summed_up, ["incoming","outgoing"])} ="incoming" 
    end
        
        % get rid of zero entries
        non_zero_rxns = find(~all(samples == 0,2));
        % subset S and sampling matrix accordingly 
        S_matrix = S_matrix(:,non_zero_rxns);
        samples = samples(non_zero_rxns,:);
        
        fluxsum=zeros(size(S_matrix,1),size(samples,2));
        for counter=1:size(samples,2)
            v=samples(:,counter); % one sample
            temp=repmat(v',size(S_matrix,1),1); %
            fluxes=S_matrix.*temp;
            
            fluxSumP=full(sum((fluxes>0).*fluxes,2));
            fluxSumN=full(sum((fluxes<0).*fluxes,2));

            if abs(fluxSumP) ~= abs(fluxSumN)
                error("Your fluxes seem to be missassigned. There was some mix up, resulting in the incoming Fluxsum not being equal to the outgoing fluxsum. ")
            end
            if flux_summed_up == "outgoing"
                fluxsum(:,counter)=fluxSumN;
            else
                fluxsum(:,counter)=fluxSumP;
            end
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