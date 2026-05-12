function [project, comparisonName] = modelsComparison(project, modelList,analysisID,referenceModel,analyses,identifier)
    % This function runs a set of analysis for the comparison of the
    % specified models.
    % A number of analysis are run: 
    % - structural analysis: based on the differential presence of
    %   metabolites, genes and reactions in the different models 
    % - functional analysis: based on the quantitative values like FVA,
    %   FBA
    % - sampling analysis: comparative analysis based on the sampling
    %   results
    % 
    % Inputs: 
    %   - project:          the object which is the output of the single_model_analysis
    %                       entailing the results of fba,fva,sampling, single gene
    %                       deletion etc. for a single model 
    %   - modelList:        the list of Model names to be included in the comparison 
    %   - analysisID:       the individual analysis slots used for the
    %                       comparison (as a default the analysis performed
    %                       most recent is choosen)
    %   - referenceModel:  the reference model used to compute the relative reaction presence
    %   - analyses:         the list of analyses which should be performed 
    %                       + modelStructureComparison: investigates the differences between the models 
    %                         on a structural level, is a gene,metabolite,or rxns present or not ? 
    %                       + modelFunctionalComparison: investigates the differences between models on a functional
    %                         level, how much flux do reactions carry in FVA, FBA solutions? 
    %                       + aamplingComparison: investigates the differences between models sampling solutions,
    %                         investigates samples solution space 
    %                         
    %   - identifier:       a string, will be added as a postfix to the analysis name, can be choosen freely
    %                       default)
    %   
    % Output : 
    %   - project:          project object with a added comparison field entailing
    %                       all the output, modelcomparison information
    %   - comparisonName:   gives back the name of the comparison name added
    arguments
        project (1,1) struct %{assertValidProjectStruct(project)}
        modelList (1,:) string
        analysisID (1,:) string
        referenceModel (1,1) string = "orig_model"
        analyses  (1,:) string  {mustBeMember(analyses, ["modelStructureComparison","modelFunctionalComparison","modelSamplingComparison","IDAREoutput"])} =["modelStructureComparison"]
        identifier (1,1) string = string(datetime('now','Format','_yyyyMMdd_HHmmss'))
    end
    
    validateInput(project, modelList, analysisID, referenceModel)

    % -- get mapping status 
    % -> mapping of rxns based on the discretization of the gene expression
    % data, as active or n
    project.models = structfun(@(x) getMappedStatus(x), ...
                          project.models, ...
                          'UniformOutput', false);
    % -- reordering discretized data matrix 
    % reorder the discretized expression data matrix in every model to
    % match the order of the genes in each model + add gene symbols
    project.models = structfun(@(x) reorderDiscretizedToMatchGeneOrder(x), ...
                              project.models, ...
                              'UniformOutput', false);

    % -- create comparison slot
    % give the comparison the name of all models + a identifier choosen

    order = string(fieldnames(project.models));
    [~, idx] = ismember(modelList, order);
    [~, sortIdx] = sort(idx);
    
    modelListOrdered = modelList(sortIdx);
    clear modelList
    comparisonName = join(modelListOrdered, "_vs_") + "__" + identifier;
    % does this comparison object already exist, was the structural
    % analysis performed ? The structural analysis is needed to do the two
    % other analysis, therefore we check if it was already run, to not run
    % it again if it was already created
    if isstruct(project.comparisons)
        fields = fieldnames(project.comparisons);
    else
        disp('project.comparisons is empty or not a struct');
        fields = {};
    end

    if ismember(comparisonName,fields)
        structureAnalysisAlreadyRun = project.comparisons.(comparisonName).structure_analysis_status;
    else
        structureAnalysisAlreadyRun = 0;
    end
    
    % run structure comparison if: it was the input of the function, if
    % another comparison was input but the structure analysis was not run
    % yet (no comparison object in the comparisons slot, or a comparison
    % object that is not complete -> structure_analysis_already_ran  not
    % defined
    if any(matches(analyses, "modelStructureComparison")) | ~structureAnalysisAlreadyRun | ~ismember(comparisonName,string(fieldnames(project.comparisons)))
        project.comparisons.(comparisonName).modelList = modelListOrdered;
        project.comparisons.(comparisonName).referenceModel = referenceModel;
        project.comparisons.(comparisonName).analysisID = analysisID;
        
        % -- run structural comparison - always has to be run 
    
        project.comparisons.(comparisonName) = modelStructuralComparison(project,modelListOrdered,referenceModel);
        project.comparisons.(comparisonName).referenceModel = referenceModel;
        project.comparisons.(comparisonName).structure_analysis_status = 1;
    end

    % -- run functional comparions
    if any(matches(analyses, "modelFunctionalComparison"))
        project.comparisons.(comparisonName).plots.funct = modelFunctionalComparison(project, comparisonName);
    end

    % -- run sampling comparison
    if any(matches(analyses, "modelSamplingComparison"))
        project = modelSamplingComparison(project,comparisonName);
    end
   


    function validateInput(project, modelList, analysisID, referenceModel)

        % This function checks that the input is valid and the comparison
        % can run through!
        % check that all the specified model names(modelList & referenceModel)
        validModels = string(fieldnames(project.models));
        modelsToTest = [modelList,referenceModel]
        if ~all(ismember(modelsToTest, validModels))
            invalid = modelsToTest(~ismember(modelsToTest, validModels));
            error("Either for the modelList or the referenceModel invalid model(s) have been choosen.Invalid model name(s): %s. Valid models are: %s", ...
                   strjoin(invalid,", "), ...
                   strjoin(validModels,", "));
        end
        % check that every model has a analysisID assigned
        assert(length(modelList)== length(analysisID), 'The number of models given to analyse does not match the number of analysisIDs given as input, every model needs to have one analysisID on basis of which the comparative analysis is performed!')
        
        % check that the assigned analysisID is also present in the
        % models.analysis slot 
        modelsCheckAnalysisID = rmfield(project.models, setdiff(fieldnames(project.models), modelList));

        analysisDone = structfun(@(m)fieldnames(m.analysis), modelsCheckAnalysisID, 'UniformOutput', false);
        for m=1:numel(modelList)
           assert(ismember(analysisID(m),analysisDone.(modelList(m))),...
                  'Model %s does not have an analysis object with the id %s',...
                  modelList(m),analysisID(m));
        end

    end

    function modelStruct = getMappedStatus(modelStruct)
        
        % This function gives you the reactions which were part of the
        % initialCore in the construction of the mdoel (1) and the rxns which
        % were discretized to be notExpressed (-1).
        % The function steps in every single model of the project object and
        % checks if there is a discretized expression data slot in the object. 
        % If the model has this data slot the the discretized data is mapped to
        % the reaction by this function and returns the mapping for every
        % colum in the discretized dataframe as well as a global mapping per
        % rxns for the whole model
    
        % - only perform the mapping for the models having associated
        % expression data 
        if any(contains(fieldnames(modelStruct),"discretized_data"))
            if isfield(modelStruct, "settings")
                if isfield(modelStruct.settings, "dico")
                    % using rfastcormics function to map discretized data to the rxns
                    mapping = mapExpressionToModel( ...
                        modelStruct.model, ...
                        modelStruct.discretized_data.values, ...
                        modelStruct.settings.dico, ...
                        string(modelStruct.discretized_data.gene_names), ...
                        1);
                    
                    numberOfSamples = size(mapping, 2);
                    % store it per sample, column in the discretized expression matrix
                    modelStruct.mappedDiscRxnsSample = mapping;
                    % and also as global mapping for the model, by multiplying it with
                    % the consensus porportion used for the model construction 
                    % parameters used in the model construction can be accessed in the 
                    % settings. slot of each individual model
                    
                    if isfield(modelStruct.settings,"script_parameters")
                        if isfield(modelStruct.settings.script_parameters,"consensus_proportion")
                            % definition of initialCore reactions
                            modelStruct.mappedDiscRxns = sum(mapping == 1, 2) >= (modelStruct.settings.script_parameters.consensus_proportion * numberOfSamples);
                            % definition of the notExpressed genes
                            notExpressed = find(sum(mapping == -1, 2) >= (modelStruct.settings.script_parameters.consensus_proportion * numberOfSamples));
                            modelStruct.mappedDiscRxns = int32(modelStruct.mappedDiscRxns);
                            modelStruct.mappedDiscRxns(notExpressed) = -1;
                            modelStruct.mappedDiscRxns = int8(modelStruct.mappedDiscRxns); % needs less storage
                        else
                            modelStruct.mappedDiscRxns = int8(repmat(-13,size(mapping,1),1));
                        end
                    else
                        modelStruct.mappedDiscRxns = int8(repmat(-13,size(mapping,1),1));
                    end
                    % The definition of the unexpressed and initialCore rxns is done as
                    % performed in rFASTCORMICS_v2
                end
            end

        end
    end

    function modelStruct = reorderDiscretizedToMatchGeneOrder(modelStruct)
        % This function reorders the discretized data slot of every model to
        % have the rows/genes in the same order as in the associated
        % model.genes slot. In addition it adds the gene symbol to the discretized data slot.
        % The discretized data is given back as a matrix. 
        if any(contains(fieldnames(modelStruct),"discretized_data"))
            if isfield(modelStruct,"settings")
                if isfield(modelStruct.settings,"dico")
                    if size(modelStruct.discretized_data,1) ~= length(string(modelStruct.settings.dico.SYMBOL))% in case this is the second time you run a modelComparison, the disc slot is already ordered so we skip this function in that case
                        geneSymbol = string(modelStruct.settings.dico.SYMBOL);
                        geneIdInModel = string(modelStruct.settings.dico.gene_id_in_model);
                        discTbl     = modelStruct.discretized_data;   % table with gene_names + data
                        % Map discretized genes into full gene list
                        [isPresent, idx] = ismember(geneSymbol, string(discTbl.gene_names));
                        % Preallocate output table with NaNs
                        outTbl = zeros(numel(geneSymbol), size(discTbl.values,2));
                        % Fill rows that exist
                        outTbl(isPresent,:) = discTbl{idx(isPresent), "values"};
                        % Add gene names as first column
                        modelStruct.discretized_data = table(geneIdInModel,geneSymbol,outTbl, 'VariableNames',...
                                                             ["gene_id_in_model",string(modelStruct.discretized_data.Properties.VariableNames)]);
                    end
                end
            end
        end
    end


end

function plots = modelFunctionalComparison(project, comparisonName)
    % This function runs the functional model comparison. 
    % The models are compared on basis of the FBA & FVA results from the singleModelAnalysis. 
    % So the functional capacity the model has in context of the defined
    % objective function. 
    % Input:
    %   - project: struct with content defined in the README,
    %     singleModelAnalysis run on the object, and chooseActiveAnalysis
    %     needs to be set too
    %   - comparisonName: which of the comparisons should be visualized,
    %     comparison slot is created when running the
    %     modelStructuralComparison function
    % Output: #TODO
    %   - different figures
    %   - excel sheet - imports export fba, rxns wise FVA scores for each model, enrichment
    %     table with the scores -> needs to be done still 
    arguments
        project (1,1) struct
        comparisonName (1,1) string
    end

    modelList = project.comparisons.(comparisonName).modelNames;
    referenceModel = project.comparisons.(comparisonName).referenceModel;
    
    
    %%% ---------- Visualization: objective values per model
    fbaObjectiveValues = cell2mat(cellfun(@(x) project.models.(x).analysis.active.FBA.f(1,1) ,modelList,"UniformOutput",false));
    getExchangeRxnsIdx = find(findExcRxns(project.models.(referenceModel).model));    

    plots.objValue = figure('Color','w','Position',[20 20 700 300],'Visible','off');
 
    bar(fbaObjectiveValues)
    title('Model comparison: flux of optimized reaction')
    ylabel('Reaction flux value for objective function [mMol/(gDW*h)]')
    xlabel('Model')
    xticklabels(modelList)
    set(gca, 'FontSize', 18)

    %%% ---------- Visualization: get FBA values under objective function
    %%%             for the different models - filtered for exchange rxns
    
    % Import
    plots.import = getFluxPlot(project, comparisonName,getExchangeRxnsIdx, ...
                                    'thresholdFlux','upper','FVA',false,'reducedCost',false,'visiblePlots',"off");
    % Export
    plots.export = getFluxPlot(project, comparisonName,getExchangeRxnsIdx,...
                  'thresholdFlux','lower','FVA',false,'reducedCost',false,'visiblePlots',"off");
    
    %%% ---------- Visualization: FVA Similarity between Models

    [fvaSim,fvaSimRxns, fvaSimPathways] = computeFvaSimilarity(project,comparisonName);

    plots.fvaSim.overall = plotClustergram(fvaSim,...
                             modelList,...
                             modelList,...
                             {'Similarity of FVA boundaries'},...
                             "FVA similarity",...
                             [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255);
    
    %%% ---------- Visualization: FVA Similarity per reaction histogramm

    plots.fvaSim.hist = FVASimValuesHist(fvaSimRxns, modelList);
    
    %%% ---------- Visualization: FVA Similarity per reaction - enrichment
    %%%            for low fva similarity scores per pathway in the model



    resEnrichment = getEnrichmentTable(project,modelList,fvaSimRxns,referenceModel,[]);
    % put the results of FDR and NES in one matrix each

    comparisons = fieldnames(resEnrichment);

    % All unique pathways
    allPathways = unique(vertcat(resEnrichment.(comparisons{1}).Subsystem));
    for k = 2:numel(comparisons)
        allPathways = unique([allPathways; resEnrichment.(comparisons{k}).Subsystem]);
    end

    % Preallocate tables
    NESTbl = array2table(nan(numel(allPathways),numel(comparisons)), ...
        'RowNames', allPathways, 'VariableNames', comparisons);
    FDRTbl = array2table(nan(numel(allPathways),numel(comparisons)), ...
        'RowNames', allPathways, 'VariableNames', comparisons);

    % Fill tables
    for c = 1:numel(comparisons)
        comp = comparisons{c};
        idx = ismember(allPathways, resEnrichment.(comp).Subsystem);
        [~, loc] = ismember(allPathways(idx), resEnrichment.(comp).Subsystem);
        NESTbl{idx, comp} = resEnrichment.(comp).NES(loc);
        FDRTbl{idx, comp} = resEnrichment.(comp).FDR(loc);
    end

    filterForSig = find((sum(FDRTbl{:,:} < 0.05,2) > 0) & (sum(NESTbl{:,:} > 0,2) >0));
    FDRTbl = FDRTbl(filterForSig,:);
    NESTbl = NESTbl(filterForSig,:);

    plots.fvaSim.enrich = dotplot(NESTbl,FDRTbl);
    
    
    %%% ---------- Visualization: Fluxsum based on the FBA values ? 

    replacementValue = "analysis.active.FBA.v"; % get the fba solution values
    project.comparisons.(comparisonName).orderedFba = getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacementValue);
    
    % compute Fluxsum 

    [idxPathways,names_pathways] = getDefaultSubsystems(project, referenceModel);                                 

    [fluxsumSets,plots.fba.heatmapRxnFluxsum] = visualizeFluxsum(project,comparisonName,[],idxPathways,...
                                                                     names_pathways,...
                                                                     "heatmap",true,true,"orderedFba", "reactions",...
                                                                      referenceModel,"off");
    
    
    plots.fba.heatmapRxnActivityFba = getNetworkActivity(project,comparisonName,idxPathways,names_pathways);


    [fluxsumSets,plots.fba.heatmapMetsFluxsum] = visualizeFluxsum(project,comparisonName,[],idxPathways,...
                                                                      names_pathways,...
                                                                      "heatmap",true,true,"orderedFba", "incoming",...
                                                                      referenceModel,"off");

    %%% -> show the top 20 most variant metabolites excluding known cofactors 
    % cofactorNames = ["atp", "adp", "amp", "nad", "nadh", "nadp", "nadph", ...
    %             "coa", "accoa", "fad", "fadh2", "pi", "pp_i"];

    %%% -> show the fluxsum in the defined pathways
    % get_essential_pathway_metabolites('Glycolysis',project,referenceModel)

    function fig = FVASimValuesHist(fvaSimRxns, modelList)
        % This function visualizes the FVA values in a histogramm per
        % comparison. These histogramms give us an indication of how similary
        % models are in their FVA min and max boundaries, and although the
        % heatmap summing up the FVA values also hast the same intuition, this
        % visualization also enables to see what similarity values are mostly
        % occuring. Is the difference in the overall similarity mainly due to a
        % few reactions having very low values, or do we see a lot of mean fva
        % similarity values per reaction ? 
        % #TODO improve function documentation!!
        fvaLower2x2 = getLowerTriangleBlock(fvaSimRxns);
        
        modelPairs = cell(numel(modelList));  % preallocate
        
        for i = 1:numel(modelList)
            for j = 1:numel(modelList)
                modelPairs{i,j} = {};
                if i ~= j
                    modelPairs{i,j} = {modelList{i}, modelList{j}};
                end
            end
        end
        modelPairs2x2 = getLowerTriangleBlock(modelPairs);
        
        [nRows, nCols] = size(fvaLower2x2);
        
        fig = figure('Color','w','Visible','off','Position', [100 100 2000*3 2000]);
        % Create tiled layout
        t = tiledlayout(fig,nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
        
        for i = 1:nRows
            for j = 1:nCols
                nexttile((i-1)*nCols + j)
        
                data = fvaLower2x2{i,j};
    
                if ~isempty(data)
                    data = data(data ~= 1);  % remove trivial values
                end
            
                if ~isempty(data)
                    histogram(data,100)
                    set(gca, 'FontSize', 18)
                else
                    axis off  % empty tile
                end
        
                % Label axes
                if ~isempty(data)
                    if i ~= 1
                        xlabel(modelPairs2x2{i,j}{1,2}, 'Interpreter','none')
                    end
                    if j == 1
                        ylabel(modelPairs2x2{i,j}{1,1}, 'Interpreter','none')
                    end
                end
                box on
            end
        end
        
        sgt = sgtitle('Histogram of FVA similarity values between models per rxns (<1)');
        sgt.FontSize = 20;   % set desired font size
    end

    function Results = getEnrichmentTable(project,modelList,fvaSimRxns,referenceModel, subSystems)
        % This function visualizes the enrichment results in a dotplot!!
        % #TODO: better documentation of the function1!!!
    
        [~,rxnMapping] = getOrderedFeatureMatrix(project,modelList,"rxns", referenceModel);
    
    
        if isempty(subSystems)
            subSystems = string(project.models.(referenceModel).model.subSystems); 
        end
        rxns = string(project.models.(referenceModel).model.rxns);        
    
        [uniqSubs, ~, idx] = unique(subSystems);
    
        subStruct = struct();
    
        for k = 1:numel(uniqSubs)
            subName = uniqSubs{k};
            fieldName = matlab.lang.makeValidName(subName);
        
            % Collect reactions belonging to this subsystem
            subStruct.(fieldName).name = subName;
            subStruct.(fieldName).rxns = rxns(idx == k);
        end
        n = length(modelList);
        [I, J] = ndgrid(1:n, 1:n);
        modelIndex = arrayfun(@(i,j) [i j], I, J, 'UniformOutput', false);
    
        fvaSim = getLowerTriangleBlock(fvaSimRxns);
        modelIndex = getLowerTriangleBlock(modelIndex);
     
        modelPairs = cell(numel(modelList));  % preallocate
        
        for i = 1:numel(modelList)
            for j = 1:numel(modelList)
                modelPairs{i,j} = {};
                if i ~= j
                    modelPairs{i,j} = {modelList{i}, modelList{j}};
                end
            end
        end
        modelPairs2x2 = getLowerTriangleBlock(modelPairs);
        
        
        Results    = struct();
        for k = 1:numel(fvaSim)
    
            x = fvaSim{k};
            y = strjoin(modelPairs2x2{k},'_');
    
            model1idx = modelIndex{k}(1);
            model2idx = modelIndex{k}(2);
    
            if isempty(x) || isempty(y)
                continue
            end
    
            rxnIdsInBothModels = find(sum(rxnMapping(:,[model1idx,model2idx]) ~= 0,2) ==2);
            % filter for the rxn similarities that are in both models 
            rxnsInBothModels = rxns(rxnIdsInBothModels);
            
            Results.(string(y)) = pathwayEnrichment(subStruct , x(rxnIdsInBothModels),rxnsInBothModels);
    
        end
    
    
        function results = pathwayEnrichment(sets, metricMatrix,featureNames)
        % This function performs pathway enrichment on the fva similarity
        % values. In the context of metabolic modelling the enrichment in this
        % function is defined as the enrichment of low fva similarity values
        % (high dissimilarity) in a specific metabolic pathway. 
        % So practically this function does a ranked based hyptothesis testing.
        % Sorting of the rxns after their similarity value (ascending sorting,
        % since we want to know where the FVA boundaries are most different)
        % and then see which of the metabolic subsystems defined in the model
        % are enriched in this sorting!
        % # TODO -> better documentation of the function + check what happens
        % with the rnxs that have the same rank, there should be a group fo
        % rxns that have the same rxn FVA similarity -> does this effect the
        % enrichment -> since the sorting with the same value is then kind of
        % arbitraty !!!! 
    
    
        [metricSorted, sortIdx] = sort(metricMatrix,'descend');
        rxnsSorted = featureNames(sortIdx);
        N = numel(featureNames);
        nPerm = 1000;   % permutations
        p = 1;          % weight exponent (0 = unweighted)
        weights = abs(metricSorted).^p;
    
    
        subNames = fieldnames(sets);
        nSets = numel(subNames);
        
        ES   = nan(nSets,1);
        NES  = nan(nSets,1);
        pval = nan(nSets,1);
        setSize = nan(nSets,1);
        
        for s = 1:nSets
        
            subField = subNames{s};
            subRxns  = sets.(subField).rxns;
            setSize(s) = numel(subRxns);
        
            % Skip very small subsystems
            if setSize(s) < 5
                continue
            end
        
            % Hits in ranked list
            hits = ismember(rxnsSorted, subRxns);
            Ns = sum(hits);
        
            % ----- observed enrichment score -----
            Phit  = weights .* hits;
            Phit  = Phit / sum(Phit);
            Pmiss = (~hits) / (N - Ns);
        
            runningSum = cumsum(Phit - Pmiss);
        
            [~, imax] = max(abs(runningSum));
            ES(s) = runningSum(imax);
        
            % ----- permutation test -----
            ESnull = zeros(nPerm,1);
        
            for k = 1:nPerm
                permHits = hits(randperm(N));
        
                Phit_p  = weights .* permHits;
                Phit_p  = Phit_p / sum(Phit_p);
                Pmiss_p = (~permHits) / (N - Ns);
        
                rs_p = cumsum(Phit_p - Pmiss_p);
                ESnull(k) = max(abs(rs_p));
            end
        
            % empirical p-value
            pval(s) = mean(abs(ESnull) >= abs(ES(s)));
        
            % normalized enrichment score
            NES(s) = ES(s) / mean(abs(ESnull));
        end
    
        qval = mafdr(pval, 'BHFDR', true);
        results = table( ...
        subNames, setSize, ES, NES, pval, qval, ...
        'VariableNames', {'Subsystem','Size','ES','NES','pValue','FDR'} );
    
        results = sortrows(results, 'NES', 'descend');
    
    
        %results = results(results.FDR < 0.05, :);
    
    
        % Choose subsystem to plot
        % subField = 'HeparanSulfateDegradation'; 
        % subRxns  = sets.(subField).rxns;
        % 
        % % Compute hits
        % hits = ismember(rxnsSorted, subRxns);
        % Ns = sum(hits);
        % 
        % % Compute running sum for enrichment
        % Phit  = weights .* hits;
        % Phit  = Phit / sum(Phit);
        % Pmiss = (~hits) / (N - Ns);
        % 
        % runningSum = cumsum(Phit - Pmiss);
        % 
        % % Plot
        % figure('Color','w','Position',[200 200 600 400])
        % plot(runningSum, 'b', 'LineWidth', 2)
        % hold on
        % yline(0,'k--','LineWidth',1)
        % 
        % % Mark hits as vertical lines along x-axis
        % hitIdx = find(hits);
        % for i = 1:numel(hitIdx)
        %     x = hitIdx(i);
        %     line([x x], [min(runningSum) 0], 'Color',[0.7 0.7 0.7],'LineStyle','-')
        % end
        % 
        % xlabel('Reactions ranked by FVA similarity')
        % ylabel('Running enrichment score')
        % title(['Enrichment of ', sets.(subField).name])
        % grid on
    
    
    end
 end   
function fig = dotplot(NESTbl,FDRTbl)
    % #TODO better documentation of the function!!
    
    % --- Sort pathways by overall NES magnitude ---
    [~, sortedIdx] = sort(sum(abs(NESTbl{:,:}),2), 'descend');
    NESTbl = NESTbl(sortedIdx,:);
    FDRTbl = FDRTbl(sortedIdx,:);
    
    % --- Handle zeros in FDR and transform ---
    lowValues = 1e-10;
    FDRTbl{:,:}(FDRTbl{:,:} == 0) = lowValues;
    FDRTbl{:,:} = -log10(FDRTbl{:,:});
    
    % --- Extract labels ---
    pathways = regexprep(string(NESTbl.Properties.RowNames), "_", " ");
    comparisons = string(NESTbl.Properties.VariableNames);
    
    % --- Extract numeric matrices ---
    NES = NESTbl{:,:};
    FDR = FDRTbl{:,:};
    
    nP = numel(pathways);
    nC = numel(comparisons);
    
    % --- Create grid for scatter ---
    [X, Y] = meshgrid(1:nC, 1:nP);
    x = X(:);
    y = Y(:);
    
    nesVals = NES(:);
    fdrVals = FDR(:);
    
    % --- Prepare FDR for coloring ---
    cVals = fdrVals;              
    cVals(cVals < -log10(0.05)) = NaN;  % values >0.05 will be grey

    % --- Dot size proportional to |NES| with enhanced visual difference ---
    scatterMin = 10;    % smallest dot area (points^2)
    scatterMax = 500;  % largest dot area (points^2)

    nes = nesVals(~isnan(cVals));
    
    absNESNorm = (abs(nes) - min(abs(nes))) / (max(abs(nes)) - min(abs(nes))); % normalize 0-1
    dotSize = scatterMin + (absNESNorm.^0.5) * (scatterMax - scatterMin);  % power 0.5 emphasizes large values
    
    
    % --- Create figure ---
    fig = figure('Color','w','Position',[100 100 1000 1000],"Visible","off");
    hold on
    
    % Scatter plot for significant FDR (≤ 0.05)
    scatter(x(~isnan(cVals)), y(~isnan(cVals)), dotSize, cVals(~isnan(cVals)), 'filled')
    hold on
    
    % --- Size legend for |NES| with min, percentiles, max ---
    sizeVals = [min(abs(nes)), ...
                 prctile(abs(nes), 25), ...
                 prctile(abs(nes), 50), ...
                 prctile(abs(nes), 75), ...
                 max(abs(nes))];
    
    sizeScaled = [min(dotSize), ...
                 prctile(dotSize, 25), ...
                 prctile(dotSize, 50), ...
                 prctile(dotSize, 75), ...
                 max(dotSize)];
   
    % Custom "legend" inside axes
    legendSizes = sizeScaled;      % sizeData
    legendLabels = string(round(sizeVals,2));
    legendX = max(x) + 1;  % x position outside plot
    legendY = y(1:length(legendSizes)) ;        % y positions
    
    for i = 1:length(legendSizes)
        scatter(legendX, legendY(i), legendSizes(i), 'k', 'filled')
        text(legendX+0.2, legendY(i), legendLabels{i}, 'FontSize', 12, 'VerticalAlignment','middle')
    end
    % --- Axes formatting ---
    xticks(1:nC)
    xticklabels(regexprep(comparisons,"_", " vs "))
    yticks(1:nP)
    yticklabels(pathways)
    
    xlabel('Model comparison')
    ylabel('Pathway')
    
    xlim([0.5, nC + 1.5])
    ylim([0.5, nP + 0.5])
    set(gca,'YDir','reverse','FontSize',18)
    title("Pathway enrichment (dot size = |NES|, color = FDR)")
    
    % --- Colorbar ---
    nColors = 256;
    cmap = [linspace(1,0,nColors)' linspace(0,0,nColors)' linspace(0,1,nColors)']; 
    colormap(cmap)        % red (low) -> blue (high)
    clim([-log10(0.05) -log10(lowValues)])       
    cb = colorbar;
    cb.Label.String = '-log10(FDR)';
    cb.FontSize = 14;
    
end

function fvaFower = getLowerTriangleBlock(fvaSimRxns)
        % This function gets the lower part of a similarity matrix. 
        % Written because we are looking at pairwise distances/similarities in
        % the figures, but to not have repetitive plots it is useful to only
        % look at one part of the matrix (eiter above or below the diagonal)
        % since the diagonal (so the comparison of the model with itself) will
        % always be 1 or zero (depending wheterh we talk about distance or
        % similarity) and therefore only one of the triangels in the matrix
        % cell is interesting. 
        % This function sets all but the lower triangle to 0 so that there are
        % no repetitive plots!
        % fvaSimRxns: n x n cell array of comparisons
        % Returns a compact lower-triangle block cell array
    
        n = size(fvaSimRxns,1);
        if n ~= size(fvaSimRxns,2)
            error('Input must be a square cell array');
        end
    
        % Count how many rows/columns to keep: remove first row if needed
        % General: keep lower-triangle below the first row
        rowsToKeep = 2:n;       % always skip the first row
        colsToKeep = 1:(n-1);   % always skip the last column
    
        % Take the lower-triangle block
        fvaFower = fvaSimRxns(rowsToKeep, colsToKeep);
end


end

function structureAnalysis = modelStructuralComparison(project, modelList,referenceModel)
    % The structure comparison is a function that compares the models
    % listed on structural differences. Structural differences in the
    % context of Fastcore can be defined as the set of reactions that are
    % kept when running fastcore. This means we check for the existence of
    % rxns, metabolites and genes in the model, and their overlap between
    % models. Reaction existence will be analysed in different sets.
    % Reaction exisitence in different metabolic subsystems/pathways,
    % existence of reaction depending on their initial classification of
    % expressed/not expressed core/not expressed. This figures are meant to
    % give an indication of which reactions were kept in the model, how
    % high the overlap of those are between the models, which of the
    % overlapping reactions are within the core or not expressed, which
    % pathways the genes in the outer and intersection are part of etc.
    % Input:
    %   - project: fastcore workflow project, with a run of 
    %              singleModelAnalysis already performed, and the active
    %              analysis already set using
    %              chooseActiveAnalysisForComparison
    %   - modelList: List of models to compare
    %   - referenceModel: for the structural comparison of wether a
    %                      reaction is existent or not, a reference model needs to be defined
    %                      therefore you need to define a reference model
    %                      which is also in the models slot of the project
    %                      object
    % Output:
    %   - structure analysis: yet to be structured properly #TODO
    arguments
        project (1,1) struct
        modelList (1,:) string
        referenceModel (1,1) string
    end

    modelsList = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    models = structfun(@(x) x.model, modelsList, 'UniformOutput',false);
    structureAnalysis.modelNames = string(fieldnames(models));


    % get model sizes - # genes,reactions and metabolites
    data = struct2array(structfun(@(x) {numel(x.rxns); numel(x.mets); numel(x.genes)}, ...
                           models, 'UniformOutput', false))';
    array2table(data, ...
                    'VariableNames', {'count_reactions','count_metabolites','count_genes'}, ...
                    'RowNames', string(fieldnames(models))')

    [rxnPresence,rxnMapping] = getOrderedFeatureMatrix(project,modelList,"rxns", referenceModel);
    structureAnalysis.rxn_mapping_table = array2table(rxnMapping,"VariableNames",modelList,"RowNames",string(project.models.(referenceModel).model.rxns));

    
    % -- Visualization: Discretization status for expression of genes in model on sample level, on model level as well as the mapping/discretization on rxn level
    % -> gives you a feeling of how many reactions in the model are from the core, how many of the rxns that were notExpressed made it in regardless etc.
    
    % check if the models has all the information needed to analyse the
    % core reactions
    
    
    if all(structfun(@(x) isfield(x,"discretized_data"),modelsList))
        % get the reaction mapping (sample and model level) as well as the discretization values for each reaction/gene in the model 
        replacementValue = "mappedDiscRxnsSample"; % get the fba solution values
        orderedMappingRxnMatrixSampleWise = int8(getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacementValue));
        replacementValue = "mappedDiscRxns"; % get the fba solution values
        orderedMappingRxnMatrix = int8(getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacementValue));
        replacementValue = "discretized_data.values"; % get the fba solution values
        orderedMappingExprDiscMatrix = int8(getOrderedFeatureMatrix(project,modelList,"genes",referenceModel,replacementValue));
        
        if all(structfun(@(x) isfield(x.settings,"script_parameters"),modelsList))
            if all(structfun(@(x) isfield(x,"sample_metadata"),modelsList))
                if all(structfun(@(x) isfield(x.settings.script_parameters,"columns_to_define_model_samples_on"),modelsList))
                    % get the names of the single samples from the metadata slot - used in the following plots
                    columnNames = struct2cell(structfun(@(x)  string(x.sample_metadata{:,1}) + "_" + ...
                                                  string(x.sample_metadata.(x.settings.script_parameters.columns_to_define_model_samples_on)),...
                                            modelsList,"UniformOutput",false));
                    columnNames = vertcat(columnNames{:});
                else
                    columnNames = string(1:size(orderedMappingRxnMatrixSampleWise,2));
                end
            else
                columnNames = string(1:size(orderedMappingRxnMatrixSampleWise,2));
            end
        else
            columnNames = string(1:size(orderedMappingRxnMatrixSampleWise,2));
        end
        
        % get the data into one object to loop over
        datasets = { orderedMappingExprDiscMatrix,orderedMappingRxnMatrixSampleWise, orderedMappingRxnMatrix};   % replace with your actual dataset variables
        datasetNames = ["ordered_mapping_rxn_matrix_sample_wise", "ordered_mapping_expr_disc_matrix", "ordered_mapping_rxn_matrix"];  % optional titles
        xlabelsPlots = ["Samples", "Samples", "Models"]; 
        xticksPlots = {columnNames, columnNames, string(fieldnames(modelsList))}; 
        ylabelsPlots = ["# genes for genes which are in the context specific models", "# reactions for reactions which are in the context specific models", "# reactions for reactions which are in the context specific models"];  
        titlePlots = ["after discretization: ", "after mapping the gpr rules: ", "after applying consensus proportion"];  
    
        
    
        % Determine all unique discretization values across datasets (excluding 13)
        % the value 13 has been set to indicate that the rxn/gene is not in the
        % model, so the discretization is not shown, in these figures only the
        % discretization is shown of the genes/rxns in the model, the figures
        % for all genes, rxns are done in the QC script when preparing the data
        % for the model creation !!
        plots.dataDiscretization = getDiscretizationHist(datasets);
    end     
    if all(structfun(@(x) isfield(x,"core_reactions"),modelsList))

        %%% Visualize the core reaction per model
        data = struct2cell(structfun(@(x) [ length(x.core_reactions) - sum(ismember(x.core_reactions, x.model.rxns)); ...
                                                sum(ismember(x.core_reactions, x.model.rxns));...
                                                length(x.model.rxns) - sum(ismember(x.core_reactions, x.model.rxns));...
                                                length(x.model.rxns)], ...
                                                modelsList, 'UniformOutput', false));
        data = [data{:}];
    
        plots.coreReactions = plotCoreInModel(data,modelsList);
    
        % -- Visualization: Looking in deeper into the core reactions, the core
        % is what is defined by the data, therefore portrays the underlying
        % biological chnages, so the question is which reactions are part of
        % the outer and intersections we saw in the previous venn/intersection
        % diagramm ? are the differences in core reactions only due to
        % exchange/import ? transporters ? This should be avoided!
    
        % create an upsetr plot for the all the inter and outersections
        % filter out the main intersection -> the one with the longest name
    
        plots.coreReactionsIntersections = getUpsetPlotCore(project,modelsList,structureAnalysis);
        % this function only works with a comparison of up to 4 models!!
        % otherwise the plots get too complex!

    end


    % -- Visualization: Get the jaccard similarity on basis of the
    % gene,metabolite and reaction presence in the corresponding models
    % How similar are my models structuraly, which models are more similar
    % to each other than others ? 
   
    


    for fieldToInvestigate = ["genes", "mets", "rxns"]
        [orderedFeature, ~] = getOrderedFeatureMatrix( ...
            project, modelList, fieldToInvestigate, referenceModel);

        % Plot Venn / Heatmap of intersections based on presence
        plots.intersections.(fieldToInvestigate) =  plotFlexibleVenn( ...
                                                                    orderedFeature, ...
                                                                    structureAnalysis.modelNames, ...
                                                                    "Structural model comparison: " + fieldToInvestigate + " presence",...
                                                                    "visiblePlot","off");

        % get the jaccard distances - based on reaction presence
        % Compute Jaccard distances

        plots.jaccardDist.(fieldToInvestigate) =  plotJaccard( ...
                                                                 orderedFeature, ...
                                                                 structureAnalysis.modelNames, ...
                                                                 "Jaccard similarity of " + fieldToInvestigate + " presence (0 or 1) between models",...
                                                                 "visiblePlot","off");





    end

    % -- Visualization: Get reaction presence for each model in comparison
    % to the defined reference model -> visualization per subsystem
    % Where does the difference I see in the jaccard plot come from ? form
    % which subsystem, which subsystem is most different in pairwise
    % comparison ? 

    plots.reactionPathwayPresence = pathwayPresenceHeat(project,referenceModel);


    

    

    structureAnalysis.plots.struct = plots;


    function fig = getDiscretizationHist(datasets)
        % This function produces a nested barplot that shows how many of
        % rxns have a discretization status of 0 1 or -1.For the different
        % models but also for the different samples in the models!
        numDatasets = length(datasets);
        allValues = [];
        for k = 1:numDatasets
            allValues = union(allValues, setdiff(unique(datasets{k}), 13));
        end
        allValues = sort(allValues);  % e.g., [-1 0 1]
        
        % Define a colormap with one color per value
        cmap = lines(length(allValues));
        
        % Create figure
        fig = figure('Color','w','Position',[100 100 2000*numDatasets 2000],...
                                           'Visible','off');
        
        for k = 1:numDatasets
            dataset = datasets{k};
            xlabelPlot = xlabelsPlots(k);
            xtickPlot = xticksPlots{k};
            ylabelPlot = ylabelsPlots(k);
            titlePlot = titlePlots(k);
            
            % counts per category per sample (make sure all_values are included)
            counts = zeros(length(allValues), size(dataset,2));
            for i = 1:length(allValues)
                counts(i,:) = sum(dataset == allValues(i), 1);
            end
        
            % stacked barplot in subplot
            ax = subplot(1, numDatasets, k);
            b = bar(ax, counts', 'stacked');
        
            % Assign consistent colors
            for i = 1:length(allValues)
                b(i).FaceColor = cmap(i,:);
            end
        
            % percentages
            tot = sum(counts,1);
            pct = 100 * counts ./ tot;
            pct = round(pct);
        
            % write percentages inside bars
            for i = 1:size(counts,2)
                y0 = 0;
                for j = 1:size(counts,1)
                    if counts(j,i) > 0
                        text(i, y0 + counts(j,i)/2, ...
                            sprintf('%g%%', pct(j,i)), ...
                            'HorizontalAlignment','center', ...
                            'VerticalAlignment','middle', ...
                            'FontSize',13,'Color','w','FontWeight','bold');
                    end
                    y0 = y0 + counts(j,i);
                end
            end
        
            % axes labels and formatting
            ax.FontSize = 14;
            xlabel(xlabelPlot,'FontSize',16)
            ylabel(ylabelPlot,'FontSize',16)
            title(titlePlot,'FontSize',18)
        
            xticks(1:length(xtickPlot))
            xticklabels(regexprep(xtickPlot, "_", "-"))
        end
        
        % Single legend for the whole figure
        lgd = legend(string(allValues), 'Location','northeastoutside');
        lgd.FontSize = 14;
        lgd.Title.String = "Discretization status";
    end
    function plt = pathwayPresenceHeat(project,referenceModel)
        % This function produces a heatmap showing the rxns presence in
        % comparison to a reference model in the color scale (0-1) and the
        % absolute number of rxns per pathway in every model as text in the
        % tiles, additionally it displays the pathway size of each pathway
        % on the left side in an additional barplot!
        pathways = string(project.models.(referenceModel).model.subSystems); % get pathways from reference model
        uniquePathways = unique(pathways); 
    
        % for every pathway get the matrix identifying the rnxs from reference
        % model in this pathway
        groups = arrayfun(@(x) find(pathways == x), uniquePathways, 'UniformOutput', false);
        numGroups = length(groups);
        G = zeros(numGroups, size(orderedFeature,1));
    
        for g = 1:numGroups
            G(g, groups{g}) = 1;
        end
    
        % get the presence per subsystem in the context specific models 
        % by mulitplying the rxns presence for each subsystem (matrix ordered features) 
        % with the matrices defining the subsystem for every rxns
        
        pathwayCounts = array2table(G * orderedFeature, ...
                                     'VariableNames', structureAnalysis.modelNames,...
                                     'RowNames',string(cellstr(uniquePathways)));
        % add reference model count to be able to make a relative abundance
        pathwayCounts.referenceModel = groupcounts(pathways);
    
        relativeCounts = array2table(pathwayCounts{:,1:end-1} ./ pathwayCounts.referenceModel, ...
                                     'VariableNames', structureAnalysis.modelNames,...
                                     'RowNames',cellstr(uniquePathways));
    
        % get the idx of the most variant pathways in terms of rxns presence
        relativeCounts.row_var = var(relativeCounts{:,:}, 0, 2);
        pathwayCounts.row_var = var(pathwayCounts{:,1:end-1},0,2);
        pathwayCounts = pathwayCounts(pathwayCounts.referenceModel < 1000,:);
        % Get indices of top n highest variance rows
    
    
        pathwayCounts = sortrows(pathwayCounts,"row_var","descend");
        pathwayCounts = pathwayCounts(find(pathwayCounts.row_var ~= 0),:);
        relativeCounts = relativeCounts(pathwayCounts.Properties.RowNames,:);
        % plot top 20 most variant pathways between the choosen models
        data = relativeCounts{:,1:end-1};
        rowNames = string(relativeCounts.Properties.RowNames);
        colNames = structureAnalysis.modelNames;
    
        %%%%%%%%%%
    
        plt = figure('Color','w','Position',[20 20 700 300],'Visible','off');
        tiledlayout(1,4, ...
            'TileSpacing','compact', ...
            'Padding','compact')
    
        % ---- Bar plot (LEFT) ----
        ax1 = nexttile(1);
        barh(pathwayCounts.referenceModel)
        title('Subsystem size in the reference model')
        xlabel('# rxns in the reference model')
    
        ax1.YTick = 1:numel(rowNames);
        ax1.YTickLabel = rowNames;
        ax1.YDir = 'reverse';
        ax1.YAxisLocation = 'right';   % ⭐ labels between plots
        ax1.TickLength = [0 0];        % removes tick marks
        %ax1.XTick = [];               % remove x ticks
        ax1.YTick = 1:numel(rowNames);
        ax1.YTickLabel = rowNames;    % keep labels
        ax1.YAxisLocation = 'right';  % labels between plots
        ax1.FontSize = 12;   % bar plot labels
        % Flip the horizontal axis
        ax1.XDir = 'reverse';
    
    
    
    
        % ---- Heatmap (RIGHT, spanning 2 tiles) ----
        % z-scaling of the data -> so that the colorod
        ax2 = nexttile(2,[1 3]);
    
        imagesc(data)
        cmap = getColorPallette();
        colorbar
        title("relative counts of subsystem rxn occurence/reference model" )    % grayscale
        ax2.XTick = 1:numel(colNames);
        ax2.XTickLabel = colNames;
    
        ax2.YTick = 1:numel(rowNames);
        ax2.YTickLabel = [];    % labels shown only once
    
        ax2.TickLength = [0 0];
    
        xlabel('Models')
        ax2.FontSize = 12;   % heatmap labels
    
        % Get size of the data
        [nRows, nCols] = size(data);
    
    
        % Loop over every cell and place the absolute number from pathway_counts
        for i = 1:nRows
            for j = 1:nCols
                % You want absolute numbers, not relative counts, so use pathway_counts
                % (or multiply relative_counts by referenceModel if needed)
                value = pathwayCounts{i,j}; % +1 because first column is referenceModel
                % Place text at the center of the tile
                text(ax2, j, i, num2str(value), ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'Color','k', ...          % black text
                    'FontSize',10)
            end
        end
        % ---- Align rows ----
        linkaxes([ax1 ax2],'y')

    end
    function fig = plotCoreInModel(data,modelsList)
        % This function visualizes how many of the core rxns defined by
        % rFastcormics make it in and how many of the core rxns are
        % overlapping between the models.=
        upperData = data([3,2],:);
        categories = fieldnames(modelsList)'; 
    
        % -- Visualization: get a sense of how many of the core reactions made
        % it into your model. In theory 100 percent of the reactions should be
        % in, but in practice this is not gona happen, adjusting the
        % thresholding/discretization in the preprocessing of the data for
        % rfastcormics could change the precentages you see in the following
        % figure
    
    
    
        fig = figure('Color','w','Visible','off','Position', [100 100 1500 1500]);
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
    
        % --- first barplot
        ax1 = nexttile(1);
        bar(upperData', 'stacked')
    
    
        % Labels
        set(gca, 'XTickLabel', categories, 'FontSize', 14)  % increase tick labels font
        xlabel('Models', 'FontSize', 14)
        ylabel('# rxns', 'FontSize', 14)
    
        % Legend
        legend({"non-core reactions","core reactions"}, 'Location','northwest', 'FontSize', 14)
    
        % Title
        title('Core and non-core reactions per model', 'FontSize', 14)
    
        % ---- Compute percentages of upper stack ----
        total = sum(upperData,1);                   % total per model
        percentUpper = 100 * upperData(2,:) ./ total;  % percentage of upper bar
    
        % ---- Add text labels on top of upper bars ----
        for i = 1:size(upperData,2)  % loop over models
            % x-position is the bar center, y-position is height of lower + upper
            x = i;
            y = upperData(1,i) + upperData(2,i)/2;  % middle of upper stack
            text(x, y, sprintf('%.1f%%', percentUpper(i)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold');
        end
    
    
    
    
         % --- second barplot
         data = data(1:2,:);  % only core vs non-core counts
        ax2 = nexttile(2);
        hb = bar(data', 'stacked');
    
        % Labels
        set(gca, 'XTickLabel', categories, 'FontSize', 14)
        xlabel('Models', 'FontSize', 14)
        ylabel('# rxns', 'FontSize', 14)
    
        % Legend
        legend({ "not included","included"}, 'Location','northwest', 'FontSize', 14)
    
        % Title
        title('Core reactions used to construct the models.', 'FontSize', 14)
    
        % ---- Compute percentages of upper stack ----
        total = sum(data,1);                   % total per model
        percentUpper = 100 * data(2,:) ./ total;  % percentage of upper bar
    
        % ---- Add text labels on top of upper bars ----
        for i = 1:size(data,2)  % loop over models
            % x-position is the bar center, y-position is height of lower + upper
            x = i;
            y = data(1,i) + data(2,i)/2;  % middle of upper stack
            text(x, y, sprintf('%.1f%%', percentUpper(i)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold');
        end
    
        % Tile 3 (THIS IS THE KEY)
        ax3 = nexttile(3,[1 2]);   % column 2, span both rows
        axis(ax3,'off')
        hold(ax3,'on')
    
    
        core_reactionsIncluded = struct2cell(structfun(@(x) x.core_reactions(find(ismember(x.core_reactions, x.model.rxns)))', ...
                                                modelsList, 'UniformOutput', false));
        core_reactionsIncluded = unique([core_reactionsIncluded{:}]);
    
        corePresence = structureAnalysis.rxn_mapping_table{core_reactionsIncluded,:} ~= 0;
        [figV,idxInterOutersections,~] = plotFlexibleVenn(corePresence, structureAnalysis.modelNames, ... 
                                                            "Structural model comparison: core rxns presence","visiblePlot","off");
    
    
        if string(class(figV)) == 'matlab.ui.Figure'
    
            % Find axes inside Venn figure
            axV = findobj(figV,'Type','axes','-not','Tag','legend');
            axV = axV(1);
    
            % Copy graphics
            copyobj(allchild(axV), ax3)
    
    
    
            % Fix geometry
            axis(ax3,'tight')
            axis(ax3,'equal')
            ax3.Clipping = 'off';
            close(figV)
        else 
            % store main figure handle
            mainFig = gcf;
    
            % create Venn/heatmap figure
            [figV,idxInterOutersections,~] = plotFlexibleVenn( ...
                corePresence, structureAnalysis.modelNames, ...
                "Structural model comparison: core rxns presence");
    
            % extract data
            X = figV.XData;
            Y = figV.YData;
            C = figV.ColorData;
    
            % close the temporary figure
            close(ancestor(figV,'figure'))
    
            % activate main figure again
            figure(mainFig)
    
            % place heatmap in tile
            nexttile(3,[1 2])
            heatmap(X,Y,C)
    
            title('Structural model comparison: core rxns presence')
    
    
        end


    end
    function plt = getUpsetPlotCore(project,modelsList,structureAnalysis)
        % This figure shows the pathway of the inter and outersections
        % between the core reactions in the model. This is a quality figure
        % to controll in which subsystems we see a difference based on the
        % defined core of the different models. We want more than just the
        % Transport or exchangers to be different between models, in this
        % figure we can control that!
        % The figure only works for comparisons of 4 or less models, since
        % the Venn diagramm calculating the intersections and outersections
        % only works for 4 or less, otherwise the histogramm would get to
        % big to visualize. 
        
        core_reactionsIncluded = struct2cell(structfun(@(x) x.core_reactions(find(ismember(x.core_reactions, x.model.rxns)))', ...
                                                modelsList, 'UniformOutput', false));
        core_reactionsIncluded = unique([core_reactionsIncluded{:}]);
        corePresence = structureAnalysis.rxn_mapping_table{core_reactionsIncluded,:} ~= 0;
        [figV,idxInterOutersections,~] = plotFlexibleVenn(corePresence, structureAnalysis.modelNames, ... 
                                                            "Structural model comparison: core rxns presence","visiblePlot","off");

        if string(class(figV)) == 'matlab.ui.Figure'

            namesIntersections = fieldnames(idxInterOutersections);
            [~,allIntersection] = max(cellfun(@(x) length(x), namesIntersections));
            idxInterOutersections = rmfield(idxInterOutersections, namesIntersections(allIntersection));
            % now get the pathway of every entry
            interOutersectionsPathways = structfun(@(x) string(project.models.(referenceModel).model.subSystems(x)),...
                                                     idxInterOutersections,'UniformOutput',false);
            C = struct2cell(interOutersectionsPathways);
            uniquePathways = unique(vertcat(C{:}));
    
    
            % Preprocess pathways: collapse transport
            S = structfun(@(x) regexprep(x,"^Transport.*","Transport"), interOutersectionsPathways, 'UniformOutput', false);
            pathwaysUnique = unique(regexprep(uniquePathways,"^Transport.*","Transport"));
    
            barNames = string(fieldnames(S));
            nBars = numel(barNames);
    
            % Build count matrix
            Y = cellfun(@(b) sum(S.(b)' == pathwaysUnique, 2)', barNames', 'UniformOutput', false);
            Y = cat(1, Y{:});
    
            % Sort bars by total counts (descending)
            [~, sortIdx] = sort(sum(Y,2), 'descend');
            Y = Y(sortIdx,:);
            barNamesSorted = barNames(sortIdx);
    
            % Plot
            plt = figure('Color','w','Position',[100 100 6000 2000], 'Visible','off');
            b = bar(Y, 'stacked');
    
            % Generate a qualitative colormap with enough colors
            numColors = size(Y,2);
            % Example 20-color qualitative palette (from ColorBrewer / Tableau)
            cmap = [ ...
                166 206 227;
                31 120 180;
                178 223 138;
                51 160 44;
                251 154 153;
                227 26 28;
                253 191 111;
                255 127 0;
                202 178 214;
                106 61 154;
                255 255 153;
                177 89 40;
                141 211 199;
                255 255 179;
                190 186 218;
                251 128 114;
                128 177 211;
                253 180 98;
                179 222 105;
                252 205 229] / 255;  % Normalize 0-1
    
            % Apply colors to each category
            for k = 1:numColors
                b(k).FaceColor = cmap(mod(k-1,size(cmap,1))+1,:);
            end
    
    
            % Labels and legend
            ax = gca;
            ax.XTickLabel = regexprep(barNamesSorted, "_", " ");
            ax.FontSize = 20;
            xlabel('Model intersections/outersections','FontSize',20)
            ylabel('# Core Reactions','FontSize',20)
            title('Count of Core reactions per pathway and intersection/outersection','FontSize',20)
    
            lgd = legend(pathwaysUnique, 'Location','northeast');
            lgd.FontSize = 20;
        else
            plt = [];
        end
    end
end

function project = modelSamplingComparison(project,comparisonName)


    listModelNames = strsplit(comparisonName, "__"); 
    listModelNames = strsplit(listModelNames(1),"_vs_");
    % give the comparison the name of all models + a identifier choosen
    referenceModel = project.comparisons.(comparisonName).referenceModel;
    
    % run structural model comparison
    replacementValue = "analysis.active.sampling.samples"; % get the fba solution values
    [project.comparisons.(comparisonName).orderedSamples,~,sampleLabels] = getOrderedFeatureMatrix(project,listModelNames,"rxns",referenceModel,replacementValue);
    project.comparisons.(comparisonName).sampleModelLabels = sampleLabels;

    objectiveID = find(ismember(project.models.(project.comparisons.(comparisonName).referenceModel).model.rxns,"biomass_reaction"));

    % Normalize: divide each column (sample) by its biomass flux value
    %project.comparisons.(comparisonName).orderedSamples = project.comparisons.(comparisonName).orderedSamples ./ project.comparisons.(comparisonName).orderedSamples(objectiveID, :);


    replacementValue = "analysis.active.FBA.v"; % get the fba solution values
    project.comparisons.(comparisonName).orderedFba = getOrderedFeatureMatrix(project,listModelNames,"rxns",referenceModel,replacementValue);

    %project.comparisons.(comparisonName).ksResults = fluxKSvsControl(project.comparisons.(comparisonName).orderedSamples, sampleLabels, "WT");
    
    %project.comparisons.(comparisonName).plots.sampling = visualizeSamplingLandscape(project,comparisonName,'visiblePlot',"off");


    [idxPathways,namesPathways] = getDefaultSubsystems(project, referenceModel); 


    [fluxsumSets,project.comparisons.(comparisonName).plots.sampling.heatmapRxnFluxsum] = visualizeFluxsum(project,comparisonName,[],idxPathways,...
                                                                          namesPathways,...
                                                                          "heatmap",true,true,"orderedSamples", "reactions",...
                                                                          referenceModel,"off");



    [fluxsumSets,project.comparisons.(comparisonName).plots.sampling.heatmapMetsFluxsum] = visualizeFluxsum(project,comparisonName,[],idxPathways,...
                                                                           namesPathways,...
                                                                           "heatmap",true,true,"orderedSamples", "incoming",...
                                                                          referenceModel,"off");

    [fluxsumSets,...
     project.comparisons.(comparisonName).plots.sampling.heatmapRxnFluxsumSamples] = visualizeFluxsum(project,comparisonName,[],idxPathways,...
                                                                     namesPathways,...
                                                                     "heatmapSample",true,true,"orderedSamples", "reactions",...
                                                                          referenceModel,"off");



    [fluxsumSets,...
     project.comparisons.(comparisonName).plots.sampling.heatmapMetsFluxsumSamples] = visualizeFluxsum(project,comparisonName,[],idxPathways,...
                                                                      namesPathways,...
                                                                      "heatmapSample",true,true,"orderedSamples", "incoming",...
                                                                          referenceModel,"off");

    % [fluxsumSets,...
    %  project.comparisons.(comparisonName).plots.sampling.heatmapRxnFluxsumSamplesAllFeatures] = visualizeFluxsum(project,comparisonName,[],idxPathways,...
    %                                                                  names_pathways,...
    %                                                                  "heatmapSampleAllFeatures",true,true,"orderedSamples", "reactions",...
    %                                                                      referenceModel,"off");
    % 
    % 
    % 
    % [fluxsumSets,...
    %  project.comparisons.(comparisonName).plots.sampling.heatmapMetsFluxsumSamplesAllFeatures] = visualizeFluxsum(project,comparisonName,[],idxPathways,...
    %                                                                   names_pathways,...
    %                                                                   "heatmapSampleAllFeatures",true,true,"orderedSamples", "incoming",...
    %                                                                      referenceModel,"off");


    [fluxsumSets,fig1] = visualizeFluxsum(project,comparisonName,[],{idxPathways{1}},...
                                              namesPathways(1),"violin",true,false,"orderedSamples",...
                                              "incoming",referenceModel,"off");
    [fluxsumSets,fig2] = visualizeFlux(project,comparisonName,[],{idxPathways{1}},...
                                         namesPathways(1), "all", "off");

    project.comparisons.(comparisonName).plots.sampling = mergeStructs( ...
                                                                        project.comparisons.(comparisonName).plots.sampling, ...
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
                [~, pval(r), ksstat(r)] = kstest2(fluxControl(r,:), fluxModel(r,:));
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



function assertValidProjectStruct(project, requiredFields, singleModelAnalysisFields)
% This is a test function checking that the input project object is in the
% right format to be used in the comparative analysis
arguments
    project
    requiredFields ={'model', 'settings'}
    singleModelAnalysisFields ={'expression_data', 'discretized_data', 'sample_metadata', 'core_reactions'}
end
    models = project.models;

    assert(isstruct(project), 'Input must be a struct.');
    assert(isstruct(models), 'models inside the project struct must also be a struct.');

    names  = fieldnames(models);
    
    referenceModelStatus = structfun(@(mod) checkModelFields(mod, requiredFields,singleModelAnalysisFields), ...
              models,'UniformOutput',true);

    assert(sum(~referenceModelStatus) > 2,...
           'You have %s models that are entailing all the fields to be able to analyse them with the workflow, check your object again and make sure that it applies with the structure of a singleModelAnalysis',...
            string(sum(~referenceModelStatus)));


    function referenceModelStatus = checkModelFields(mod, requiredFields,singleModelAnalysisFields)

        f = fieldnames(mod);
        missing = requiredFields(~ismember(requiredFields, f));
    
        assert(isempty(missing), ...
            'One Model is missing required fields: %s', ...
             strjoin(missing, ', '));

        
        assert(ismember("dico",fieldnames(mod.settings)), ...
               'Dico is missing from the model.settings struct, the dico is needed to map the gene IDs from the gpr rules to symbols!')

        missing = singleModelAnalysisFields(~ismember(singleModelAnalysisFields, fieldnames(mod)));

        if length(missing) ~= length(singleModelAnalysisFields)
            assert(isempty(missing), ...
                   'One Model is missing required fields: %s', ...
                   strjoin(missing, ', '));
            referenceModelStatus = 0;
        else
            %disp("This model does not have any of fields which are defined for a project, therefore we assume this is a reference model!")
            referenceModelStatus = 1;
        end
    end
end





