function [fvaSimilarity, fvaSimilarityRxns, fvaSimilarityPathways] = computeFvaSimilarity(project, comparison)
            
            modelList = project.comparisons.(comparison).modelNames;
            referenceModel = project.comparisons.(comparison).referenceModel;

            replacementValue = "analysis.active.FVA.minMaxFluxes.maxFlux";
            orderedFvaMaxMatrix = getOrderedFeatureMatrix(project, modelList, referenceModel, "rxns", replacementValue);
            replacementValue = "analysis.active.FVA.minMaxFluxes.minFlux";
            [orderedFvaMinMatrix, rxnMapping] = getOrderedFeatureMatrix(project, modelList, referenceModel, "rxns", replacementValue);
            
            n = numel(modelList);

            fvaSimilarity = eye(n); % Diagonal is 1
            fvaSimilarityRxns = cell(n);
            
            fvaData = cell(1, n);
            for i = 1:n
                fvaData{i} = [orderedFvaMinMatrix(:, i), orderedFvaMaxMatrix(:, i)];
            end

            indexPairs = find(triu(true(n), 1)); % Upper triangle linear indices
            [rowIdx, colIdx] = ind2sub([n, n], indexPairs);

            for k = 1:length(indexPairs)
                y = rowIdx(k);
                x = colIdx(k);

                [~, rxnSim] = FVAsimilarity(fvaData{y}, fvaData{x});

                % filter out reactions that are not in both of the models
                % the differences on structural level are not interesting
                % in this figure
                rxnIdxInBothModels = find(sum(rxnMapping(:, [x,y]) ~= 0, 2) == 2);
                overallSim = mean(rxnSim(rxnIdxInBothModels, :));

                % Fill both [y, x] and [x, y] due to symmetry
                fvaSimilarity(y, x) = overallSim;
                fvaSimilarity(x, y) = overallSim;

                fvaSimilarityRxns{y, x} = rxnSim;
                fvaSimilarityRxns{x, y} = rxnSim;
            end

            nModels = size(fvaSimilarityRxns, 2);

            fvaSimilarityPathways = cell(nModels);
            
            pathways = string(project.models.(referenceModel).model.subSystems); % get pathways from reference model
            uniquePathways = unique(pathways); 
        
            % for every pathway get the matrix identifying the rnxs from reference
            % model in this pathway
            groups = arrayfun(@(x) find(pathways == x), uniquePathways, 'UniformOutput', false);
            numGroups = length(groups);
        
            indexPairs = find(triu(true(nModels), 1)); % Upper triangle linear indices
            [rowIdx, colIdx] = ind2sub([nModels, nModels], indexPairs);
        
            for k = 1:length(indexPairs)
                        y = rowIdx(k);
                        x = colIdx(k);
        
                        matrixFvaRxns = fvaSimilarityRxns{y, x};
                        G = {};
                    
                        for g = 1:numGroups
                            G{g} = matrixFvaRxns(intersect(rxnIdxInBothModels, groups{g}), 1);
                        end
                        fvaSimilarityPathways{y, x} = cellfun(@mean, G);
                        fvaSimilarityPathways{x, y} = cellfun(@mean, G);
            end
end

