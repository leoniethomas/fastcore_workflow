function project = samplingComparison(project, comparisonName)


    listModelNames = strsplit(comparisonName, "__"); 
    listModelNames = strsplit(listModelNames(1), "_vs_");
    % give the comparison the name of all models + a identifier choosen
    referenceModel = project.comparisons.(comparisonName).referenceModel;
    
    % run structural model comparison
    replacementValue = "analysis.active.sampling.samples"; % get the fba solution values
    [project.comparisons.(comparisonName).orderedSamples, ~, sampleLabels] = getOrderedFeatureMatrix(project, listModelNames, "rxns", referenceModel, replacementValue);
    project.comparisons.(comparisonName).sampleModelLabels = sampleLabels;

    % Normalize: divide each column (sample) by its biomass flux value
    % objectiveID = find(ismember(project.models.(project.comparisons.(comparisonName).referenceModel).model.rxns, "biomass_reaction"));
    % project.comparisons.(comparisonName).orderedSamples = project.comparisons.(comparisonName).orderedSamples ./ project.comparisons.(comparisonName).orderedSamples(objectiveID, :);


    replacementValue = "analysis.active.FBA.v"; % get the fba solution values
    project.comparisons.(comparisonName).orderedFba = getOrderedFeatureMatrix(project, listModelNames, "rxns", referenceModel, replacementValue);


    [idxPathways, namesPathways] = getDefaultSubsystems(project, referenceModel); 


    [fluxsumSets, project.comparisons.(comparisonName).plots.sampling.heatmapRxnFluxsum] = visualizeFluxsum(project, comparisonName, [], idxPathways,...
                                                                          namesPathways,...
                                                                          "heatmap", true, true, "orderedSamples", "reactions",...
                                                                          referenceModel, "off");



    [fluxsumSets, project.comparisons.(comparisonName).plots.sampling.heatmapMetsFluxsum] = visualizeFluxsum(project,comparisonName, [], idxPathways,...
                                                                           namesPathways,...
                                                                           "heatmap", true, true, "orderedSamples", "incoming",...
                                                                           referenceModel, "off");

    [fluxsumSets, ...
     project.comparisons.(comparisonName).plots.sampling.heatmapRxnFluxsumSamples] = visualizeFluxsum(project, comparisonName, [], idxPathways,...
                                                                     namesPathways,...
                                                                     "heatmapSample", true, true, "orderedSamples", "reactions",...
                                                                     referenceModel, "off");



    [fluxsumSets, ...
     project.comparisons.(comparisonName).plots.sampling.heatmapMetsFluxsumSamples] = visualizeFluxsum(project, comparisonName, [], idxPathways,...
                                                                      namesPathways,...
                                                                      "heatmapSample", true, true, "orderedSamples", "incoming",...
                                                                      referenceModel, "off");


    [fluxsumSets, fig1] = visualizeFluxsum(project, comparisonName, [], {idxPathways{1}},...
                                              namesPathways(1), "violin", true, false, "orderedSamples",...
                                              "incoming", referenceModel, "off");
    [fluxsumSets, fig2] = visualizeFlux(project,comparisonName, [], {idxPathways{1}},...
                                         namesPathways(1), "all", "off");

    project.comparisons.(comparisonName).plots.sampling = mergeStructs(project.comparisons.(comparisonName).plots.sampling, ...
                                                                        fig1, fig2);

    
    function out = mergeStructs(varargin)
        out = struct();
        for k = 1:nargin
            f = fieldnames(varargin{k});
            for i = 1:numel(f)
                out.(f{i}) = varargin{k}.(f{i});
            end
        end
    end

    function results = fluxKSvsControl(fluxMatrix, sampleLabels, controlLabel)
        % FLUXKSVSCONTROL  KS test of each model vs control model, per reaction
        %
        % INPUTS:
        %   fluxMatrix   - nRxns x nSamples matrix (CHRR samples in columns)
        %   sampleLabels - 1 x nSamples cell array of model names per sample
        %                  e.g. {'modelA','modelA','modelB','modelB','control',...}
        %   controlLabel - string specifying which label is the control
        %                  e.g. 'control'
        %
        % OUTPUT:
        %   results - struct array, one entry per non-control model, each with:
        %             .model      model name
        %             .ksstat     nRxns x 1 KS statistic (effect size, 0-1)
        %             .pval       nRxns x 1 raw p-values
        %             .pval_adj   nRxns x 1 BH-adjusted p-values
        %             .meanDiff   nRxns x 1 mean(model) - mean(control)
        %             .signedKS   nRxns x 1 signed KS stat (direction + magnitude)
        
        
        % --- Control samples ---
        ctrlIdx     = strcmp(sampleLabels, controlLabel);
        fluxControl = fluxMatrix(:, ctrlIdx);
        
        % --- Non-control models ---
        modelNames = unique(sampleLabels(~ctrlIdx));
        nRxns      = size(fluxMatrix, 1);
        nModels    = numel(modelNames);
        
        % --- Preallocate output matrices ---
        ksMatrix       = zeros(nRxns, nModels);
        pvalMatrix     = zeros(nRxns, nModels);
        signedKSMatrix = zeros(nRxns, nModels);
        
        for m = 1:nModels
            fluxModel = fluxMatrix(:, strcmp(sampleLabels, modelNames{m}));
        
            ksstat = zeros(nRxns, 1);
            pval   = zeros(nRxns, 1);
        
            for r = 1:nRxns
                [~, pval(r), ksstat(r)] = kstest2(fluxControl(r, :), fluxModel(r, :));
            end
        
            meanDiff = mean(fluxModel, 2) - mean(fluxControl, 2);
        
            ksMatrix(:, m)       = ksstat;
            pvalMatrix(:, m)     = mafdr(pval, 'BHFDR', true);
            signedKSMatrix(:, m) = ksstat .* sign(meanDiff);
        end

        results.ksMatrix       = ksMatrix;
        results.pvalMatrix     = pvalMatrix;
        results.signedKSMatrix = signedKSMatrix;
        results.modelNames     = modelNames;
    end
end