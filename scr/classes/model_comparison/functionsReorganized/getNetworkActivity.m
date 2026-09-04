function fig = getNetworkActivity(project, comparisonName, idxPathways, namesPathways)
    arguments
        project 
        comparisonName 
        idxPathways 
        namesPathways 
    end
    
    modelNames = project.comparisons.(comparisonName).modelNames;
    mappingTable = project.comparisons.(comparisonName).structuralComparison.rxnMappingTable{:, :} ~= 0;
    orderedFba = findFieldRecursive(project.comparisons.(comparisonName), "orderedFba");

    result = reshape(cell2mat(cellfun(@(x) ...
                                      round( ...
                                            sum(orderedFba(x,:) ~= 0, 1) ./ ...
                                            sum(mappingTable(x,:), 1) * 100 ...
                                       ), ...
                                       idxPathways, ...
                                       'UniformOutput', false)...
                               )...
                      , size(mappingTable,2), [])';

    fig = figure('Color', 'w', 'Position', [100 100 800 800], 'Visible', 'off');
    imagesc(result)
    
    cmap = getColorPallette();
    h = colorbar; 
    clim([0, 100])
    
    ylabel(h, 'Percentage of active network rxns under FBA solution', 'FontSize', 18)        % Set title/label of colorbar
    axis equal tight            % Make cells square and remove extra space
    set(gca, 'XTickLabel', modelNames, 'XTick', 1:length(modelNames), ...
         'YTick', 1:length(namesPathways), 'YTickLabel', namesPathways)
    xtickangle(45)
    ax = gca;
    ax.FontSize = 18;  
    xlabel('Model', 'FontSize', 18)       
    ylabel('Reaction set', 'FontSize', 18)    
    [nRows, nCols] = size(result);
    
        % Loop over every cell and place the absolute number from pathway_counts
        for i = 1:nRows
            for j = 1:nCols
                value =  result(i, j); % +1 because first column is reference_model
                % Place text at the center of the tile
                text(j, i, num2str(value), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Color', 'k', ...          % black text
                    'FontSize', 18)
            end
        end

end