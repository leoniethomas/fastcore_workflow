function [kldOrderedResults] = performKLDivergenceComparison(kldOrderedResults,numParallelWorkers)
        % This function is meant to check if the kldivergence is bigger
        % between samples from different models than between the samples
        % from the same model. For more detailed description see
        % function: performKLDivergenceAnalysis.m
        % Input: 
        %       - kldOrderedResults     struct with subfields 
        %                - orderedkldSets:      [#rxns in reference × #of samples from all sets for all models]
        %                - modelLabels:         (1×#of samples from all sets for all models string) samples of each set is from 
        %                - interModelKld:       [#rxns in reference model×#of between the sets of one model
        %                - interModelLabels:    labels for the interModelKld (1×# of sets in all models string)
        %                - orderedkldSetLabels: label for each sample from which set it is coming from (1×#number of samples of all sets and models taken togther string)
        % 
        %       - numParrallelWorkers   numbers of parallel workers you
        %                               would like to use (3 workers should 
        %                               be fine for every laptop).
        % Output: 
        %       - kldOrderedResults:    adding the pvalue for every rxns
        %                               given an estimate if the difference between the sets from
        %                               different models is significantly bigger than the
        %                               difference of sets within the same model. 
        %                               [#rxns in reference model x #comparisons]
        %                               + adding the labels for which to
        %                               know which column in the intraModelKLD is coming 
        %                               from which comparison.
        %                   
        %                               
        arguments
            kldOrderedResults 
            numParallelWorkers =10
        end

        modelList = unique(kldOrderedResults.modelLabels);
        nModels   = numel(modelList);
        pairs     = nchoosek(1:nModels, 2);  % all unique pairs indices

        kldOrderedResults.intraModelKld = [];
        kldOrderedResults.intraModelKldLabels = [];
        
        for p = pairs'
        
            model1 = modelList(p(1));
            model2 = modelList(p(2));
            
            uniqueSetnames1 = unique(kldOrderedResults.orderedkldSetLabels(kldOrderedResults.modelLabels == model1));
            uniqueSetnames2 = unique(kldOrderedResults.orderedkldSetLabels(kldOrderedResults.modelLabels == model2));
            setLabels = [uniqueSetnames1,uniqueSetnames1];
            modelLabels = [repmat(model1,1,length(uniqueSetnames1)),repmat(model2,1,length(uniqueSetnames2))];
    
            A = 1:length(uniqueSetnames1);
            B = length(uniqueSetnames1)+1:length(uniqueSetnames2) + length(uniqueSetnames1);
            [Agrid, Bgrid] = ndgrid(A, B);
            pairs = [Agrid(:), Bgrid(:)];
            pairCell = num2cell(pairs, 2); % use as an input for arrayfun to run over all possible pairs between the n sets choosen
            pairwiseKdl = cell(length(pairs),1);
    
            delete(gcp('nocreate'))
            parpool(numParallelWorkers);
            parfor x=1:numel(pairCell)  
                setFromModel1 = pairCell{x}(1);
                setFromModel2 = pairCell{x}(2);
    
                set1Idx = kldOrderedResults.orderedkldSetLabels == setLabels(setFromModel1) & ...
                            kldOrderedResults.modelLabels == modelLabels(setFromModel1);
                set2Idx = kldOrderedResults.orderedkldSetLabels == setLabels(setFromModel2) & ...
                            kldOrderedResults.modelLabels == modelLabels(setFromModel2);
                
                pairwiseKdl{x} = getKldValuePairs(kldOrderedResults.orderedkldSets(:,set1Idx),kldOrderedResults.orderedkldSets(:,set2Idx));                          
            end
            pairwiseKdl = cell2mat(pairwiseKdl)'; 
            kldWithinModel1 = kldOrderedResults.interModelKld(:,kldOrderedResults.interModelLabels == model1);
            kldWithinModel2 = kldOrderedResults.interModelKld(:,kldOrderedResults.interModelLabels == model2);
    
            p_value_kdl1 = arrayfun(@(i) localRanksum(pairwiseKdl(i,:), kldWithinModel1(i,:)), ...
                             1:size(pairwiseKdl, 1))';
            p_value_kdl2 = arrayfun(@(i) localRanksum(pairwiseKdl(i,:), kldWithinModel2(i,:)), ...
                             1:size(pairwiseKdl, 1))';
           
            %p_adj_kdl = mafdr(p_value_kdl,'BHFDR', true);
            
            kldOrderedResults.intraModelKld = [kldOrderedResults.intraModelKld,max([p_value_kdl1,p_value_kdl2]')'];
            kldOrderedResults.intraModelKldLabels = [kldOrderedResults.intraModelKldLabels,model1 + "_vs_" + model2];
            % figure
            % histogram(p_value_kdl1, 100,'EdgeColor','b')
            % figure
            % histogram(p_value_kdl2, 100,'EdgeColor','r')
            % sum( p_value_kdl < 0.05)
            % fdr = sum( p_value_kdl < 0.05)/size(pairwiseKdl,1)
            % if fdr > 0.05
            %     error("FDR is %.1f%% — above the 5%% threshold. Sampling may not have converged.", fdr*100)
            % end
            
        end

        function p = localRanksum(a, b)
            if all(isnan(a)) || all(isnan(b))
                p = 1;
            else
                p = ranksum(a, b);
            end
        end


end

