function fig = get_network_activity(project,comparison_name,idx_pathways,names_pathways)
    arguments
        project 
        comparison_name 
        idx_pathways 
        names_pathways 
    end
    
    model_names = project.comparisons.(comparison_name).modelNames;
    mapping_table = project.comparisons.(comparison_name).rxn_mapping_table{:,:} ~=0;

    result = reshape(cell2mat(cellfun(@(x) ...
                                      round( ...
                                            sum(project.comparisons.(comparison_name).ordered_fba(x,:) ~= 0, 1) ./ ...
                                            sum(mapping_table(x,:), 1) * 100 ...
                                       ), ...
                                       idx_pathways, ...
                                       'UniformOutput', false)...
                               )...
                      , 3, [])';

    fig = figure('Color','w','Position',[100 100 800 800], 'Visible','off');
    imagesc(result)
    
    cmap = get_color_pallette();
    h = colorbar; 
    clim([0,100])
    
    ylabel(h, 'Percentage of active network rxns under FBA solution', 'FontSize', 18)        % Set title/label of colorbar
    axis equal tight            % Make cells square and remove extra space
    set(gca,'XTickLabel',model_names,'XTick',1:length(model_names), ...
         'YTick', 1:length(names_pathways), 'YTickLabel',names_pathways)
    xtickangle(45)
    ax = gca;
    ax.FontSize = 18;  
    xlabel('Model', 'FontSize', 18)       
    ylabel('Reaction set', 'FontSize', 18)    
    [nRows, nCols] = size(result);
    
        % Loop over every cell and place the absolute number from pathway_counts
        for i = 1:nRows
            for j = 1:nCols
                value =  result(i,j); % +1 because first column is reference_model
                % Place text at the center of the tile
                text(j, i, num2str(value), ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'Color','k', ...          % black text
                    'FontSize',18)
            end
        end

end