function [fva_similarity,fva_similarity_rxns,fva_similarity_pathways] = compute_fva_similariy(project,comparison)
            
            modelList = project.comparisons.(comparison).modelNames;
            reference_model = project.comparisons.(comparison).reference_model;

            replacement_value = "analysis.FVA.maxFlux"; % get the fba solution values
            ordered_fva_max_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
            replacement_value = "analysis.FVA.minFlux"; % get the fba solution values
            ordered_fva_min_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
            
            n = numel(modelList);

            
            
            fva_similarity = eye(n); % Diagonal is 1
            fva_similarity_rxns = cell(n);
            
            fvaData = cell(1, n);
            for i = 1:n
                fvaData{i} = [ordered_fva_min_matrix(:, i), ordered_fva_max_matrix(:, i)];
            end

            
            indexPairs = find(triu(true(n), 1)); % Upper triangle linear indices
            [rowIdx, colIdx] = ind2sub([n, n], indexPairs);

            for k = 1:length(indexPairs)
                y = rowIdx(k);
                x = colIdx(k);

                [overallSim, rxnSim] = FVAsimilarity(fvaData{y}, fvaData{x});
                

                % Fill both [y,x] and [x,y] due to symmetry
                fva_similarity(y,x) = overallSim;
                fva_similarity(x,y) = overallSim;

                fva_similarity_rxns{y,x} = rxnSim;
                fva_similarity_rxns{x,y} = rxnSim;
            end

            n_models = size(fva_similarity_rxns,2);

            fva_similarity_pathways = cell(n_models);
        
            
            pathways = string(project.models.(reference_model).model.subSystems); % get pathways from reference model
            unique_pathways = unique(pathways); 
        
            % for every pathway get the matrix identifying the rnxs from reference
            % model in this pathway
            groups = arrayfun(@(x) find(pathways == x), unique_pathways, 'UniformOutput', false);
            num_groups = length(groups);
        
            indexPairs = find(triu(true(n_models), 1)); % Upper triangle linear indices
            [rowIdx, colIdx] = ind2sub([n_models, n_models], indexPairs);
        
            for k = 1:length(indexPairs)
                        y = rowIdx(k);
                        x = colIdx(k);
        
                        matrix_fva_rxns = fva_similarity_rxns{y,x};
                        G = zeros(num_groups, size(matrix_fva_rxns,1));
                    
                        for g = 1:num_groups
                            G(g, groups{g}) = matrix_fva_rxns(groups{g},1);
                        end
                        fva_similarity_pathways{y,x} = mean(G,2);
                        fva_similarity_pathways{x,y} = mean(G,2);
            end
end

