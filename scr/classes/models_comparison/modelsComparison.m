function [project, comparison_name] = modelsComparison(project, modelList,analysisID,reference_model,analyses,identifier)
    % This function runs a set of analysis for the comparison of the
    % specified models.
    % A number of analysis are run: 
    % - structural analysis: based on the differential presence of
    %   metabolites, genes and reactions in the different models 
    % - functional analysis: based on the quantitative values like FVA,
    %   FBA
    % - analysis of sampling results
    % 
    % Inputs: 
    %   - project:          the object which is the output of the single_model_analysis
    %                       entailing the results of fba,fva,sampling, single gene
    %                       deletion etc. for a single model 
    %   - modelList:        the list of Model names to be included in the comparison 
    %   - analysisID:       the individual analysis slots used for the
    %                       comparison
    %   - reference_model:  the reference model used to compute the relative reaction presence
    %   - analyses:         the list of analyses which should be performed 
    %                       + modelStructureComparison: investigates the differences between the models 
    %                         on a structural level, is a gene,metabolite,or rxns present or not ? 
    %                       + modelFunctionalComparison: investigates the differences between models on a functional
    %                         level, how much flux do reactions carry in FVA, FBA solutions? 
    %                       + modelsComparisonSampling: investigates the differences between models sampling solutions,
    %                         investigates samples solution space 
    %                         
    %   - identifier:       a string, will be added as a postfix to the analysis name, can be choosen freely
    %                       default)
    %   
    % Output : 
    %   - project:          project object with a added comparison field entailing
    %                       all the output, modelcomparison information
    arguments
        project (1,1) struct
        modelList (1,:) string
        analysisID (1,:) string
        reference_model (1,1) string = "orig_model"
        analyses  (1,:) string  {mustBeMember(analyses, ["modelStructureComparison","modelFunctionalComparison","modelsComparisonSampling","IDAREoutput"])} =["modelStructureComparison"]
        identifier (1,1) string = string(datetime('now','Format','_yyyyMMdd_HHmmss'))
    end

    % --check input paramteres
    % check that all the specified model names(modelList & reference_model)
    validModels = string(fieldnames(project.models));
    modelsToTest = [modelList,reference_model]
    if ~all(ismember(modelsToTest, validModels))
        invalid = modelsToTest(~ismember(modelsToTest, validModels));
        error("Either for the modelList or the reference_model invalid model(s) have been choosen.Invalid model name(s): %s. Valid models are: %s", ...
               strjoin(invalid,", "), ...
               strjoin(validModels,", "));
    end
    
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
    comparison_name = join(modelList, "_vs_") + "__" + identifier;
    project.comparisons.(comparison_name).modelList = modelList;
    project.comparisons.(comparison_name).reference_model = reference_model;
    project.comparisons.(comparison_name).analysisID = analysisID;
    
    % -- run structural comparison - always has to be run 

    project.comparisons.(comparison_name) = modelStructuralComparison(project,modelList,reference_model);
    project.comparisons.(comparison_name).reference_model = reference_model;

    % -- fun functional comparions
    if any(matches(analyses, "modelFunctionalComparison"))
        modelFunctionalComparison(project, comparison_name);
    end

    % -- run sampling comparison
    if any(matches(analyses, "modelsComparisonSampling"))
        project = modelsComparisonSampling(project,comparison_name);
    end
    
    % -- generate output to visualize in IDARE
    if any(matches(analyses, "IDAREoutput"))
        folder_path = "./idare/";
        mkdir(folder_path);
        prepareDataForIDAREVisualization(project, comparison_name,folder_path);
    end
end

function modelFunctionalComparison(project, comparison_name)
    % This function runs the functional model comparison. 
    % The models are compared on basis of the FBA & FVA results from the singleModelAnalysis. 
    % So the functional capacity the model has in context of the defined
    % objective function. 
    % Input:
    %   - project: struct with content defined in the README,
    %     singleModelAnalysis run on the object, and chooseActiveAnalysis
    %     needs to be set too
    %   - comparison_name: which of the comparisons should be visualized,
    %     comparison slot is created when running the
    %     modelStructuralComparison function
    % Output: #TODO
    %   - different figures
    %   - excel sheet - imports export fba, rxns wise FVA scores for each model, enrichment
    %     table with the scores -> needs to be done still 
    arguments
        project (1,1) struct
        comparison_name (1,1) string
    end

    modelList = project.comparisons.(comparison_name).modelNames;
    reference_model = project.comparisons.(comparison_name).reference_model;
    
    
    %%% ---------- Visualization: objective values per model
    fba_objective_values = cell2mat(cellfun(@(x) project.models.(x).analysis.FBA.f(1,1) ,modelList,"UniformOutput",false));
    get_exchange_rxns_idx = find(findExcRxns(project.models.(reference_model).model));    
    figure
    bar(fba_objective_values)
    title('Model comparison: flux of optimized reaction')
    ylabel('Reaction flux value for objective function [mMol/(gDW*h)]')
    xlabel('Model')
    xticklabels(modelList)
    set(gca, 'FontSize', 18)

    %%% ---------- Visualization: get FBA values under objective function
    %%%             for the different models - filtered for exchange rxns
    
    % Import
    get_flux_plot(project, comparison_name,get_exchange_rxns_idx, ...
                  'threshold_flux','upper','FVA',false,'reducedCost',false);
    % Export
    get_flux_plot(project, comparison_name,get_exchange_rxns_idx,...
                  'threshold_flux','lower','FVA',false,'reducedCost',false);
    
    %%% ---------- Visualization: FVA Similarity between Models

    [fva_sim,fva_sim_rxns, fva_sim_pathways] = compute_fva_similariy(project,comparison_name);

    plot_clustergram(fva_sim,...
                             modelList,...
                             modelList,...
                             {'Similarity of FVA boundaries'},...
                             "FVA similarity",...
                             [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255);
    
    %%% ---------- Visualization: FVA Similarity per reaction histogramm

    FVA_sim_values_hist(fva_sim_rxns, modelList)
    
    %%% ---------- Visualization: FVA Similarity per reaction - enrichment
    %%%            for low fva similarity scores per pathway in the model

    [~,rxn_mapping] = getOrderedFeatureMatrix(project,modelList,"rxns", reference_model);
    idx_to_keep = find(sum(rxn_mapping ~=0,2) > 0);
    

    res_enrichment = get_enrichment_table(project,modelList,fva_sim_rxns,idx_to_keep,reference_model,[]);
    % put the results of FDR and NES in one matrix each

    comparisons = fieldnames(res_enrichment);
    
    % All unique pathways
    allPathways = unique(vertcat(res_enrichment.(comparisons{1}).Subsystem));
    for k = 2:numel(comparisons)
        allPathways = unique([allPathways; res_enrichment.(comparisons{k}).Subsystem]);
    end
    
    % Preallocate tables
    NES_tbl = array2table(nan(numel(allPathways),numel(comparisons)), ...
        'RowNames', allPathways, 'VariableNames', comparisons);
    FDR_tbl = array2table(nan(numel(allPathways),numel(comparisons)), ...
        'RowNames', allPathways, 'VariableNames', comparisons);
    
    % Fill tables
    for c = 1:numel(comparisons)
        comp = comparisons{c};
        idx = ismember(allPathways, res_enrichment.(comp).Subsystem);
        [~, loc] = ismember(allPathways(idx), res_enrichment.(comp).Subsystem);
        NES_tbl{idx, comp} = res_enrichment.(comp).NES(loc);
        FDR_tbl{idx, comp} = res_enrichment.(comp).FDR(loc);
    end

    filter_for_sig = find((sum(FDR_tbl{:,:} < 0.05,2) > 0) & (sum(NES_tbl{:,:} > 0,2) >0));
    FDR_tbl = FDR_tbl(filter_for_sig,:);
    NES_tbl = NES_tbl(filter_for_sig,:);

    dotplot(NES_tbl,FDR_tbl)

    
    %%% ---------- Visualization: Fluxsum based on the FBA values ? 
    
    % compute Fluxsum 

    %%% -> show the top 20 most variant metabolites excluding known cofactors 
    % cofactorNames = ["atp", "adp", "amp", "nad", "nadh", "nadp", "nadph", ...
    %             "coa", "accoa", "fad", "fadh2", "pi", "pp_i"];

    %%% -> show the fluxsum in the defined pathways
    % get_essential_pathway_metabolites('Glycolysis',project,reference_model)

end


function structure_analysis = modelStructuralComparison(project, modelList,reference_model)
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
    %   - reference_model: for the structural comparison of wether a
    %                      reaction is existent or not, a reference model needs to be defined
    %                      therefore you need to define a reference model
    %                      which is also in the models slot of the project
    %                      object
    % Output:
    %   - structure analysis: yet to be structured properly #TODO
    arguments
        project (1,1) struct
        modelList (1,:) string
        reference_model (1,1) string
    end

    models_list = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    models = structfun(@(x) x.model, models_list, 'UniformOutput',false);
    structure_analysis.modelNames = string(fieldnames(models));


    % get model sizes - # genes,reactions and metabolites
    data = struct2array(structfun(@(x) {numel(x.rxns); numel(x.mets); numel(x.genes)}, ...
                           models, 'UniformOutput', false))';
    array2table(data, ...
                    'VariableNames', {'count_reactions','count_metabolites','count_genes'}, ...
                    'RowNames', string(fieldnames(models))')
    
    % -- Visualization: Discretization status for expression of genes in model on sample level, on model level as well as the mapping/discretization on rxn level
    % -> gives you a feeling of how many reactions in the model are from the core, how many of the rxns that were notExpressed made it in regardless etc.
    
    % get the reaction mapping (sample and model level) as well as the discretization values for each reaction/gene in the model 
    replacement_value = "mappedDiscRxns_sample"; % get the fba solution values
    ordered_mapping_rxn_matrix_sample_wise = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
    replacement_value = "mappedDiscRxns"; % get the fba solution values
    ordered_mapping_rxn_matrix = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
    replacement_value = "discretized_data.values"; % get the fba solution values
    ordered_mapping_expr_disc_matrix = getOrderedFeatureMatrix(project,modelList,"genes",reference_model,replacement_value);

    % get the names of the single samples from the metadata slot - used in the following plots
    columnnames = struct2cell(structfun(@(x)  string(x.sample_metadata{:,1}) + "_" + ...
                                  string(x.sample_metadata.(x.settings.script_parameters.columns_to_define_model_samples_on)),...
                            models_list,"UniformOutput",false));
    columnnames = vertcat(columnnames{:});

    % get the data into one object to loop over
    datasets = { ordered_mapping_expr_disc_matrix,ordered_mapping_rxn_matrix_sample_wise, ordered_mapping_rxn_matrix};   % replace with your actual dataset variables
    dataset_names = ["ordered_mapping_rxn_matrix_sample_wise", "ordered_mapping_expr_disc_matrix", "ordered_mapping_rxn_matrix"];  % optional titles
    xlabels_plots = ["Samples", "Samples", "Models"]; 
    xticks_plots = {columnnames, columnnames, string(fieldnames(models_list))}; 
    ylabels_plots = ["# genes for genes which are in the context specific models", "# reactions for reactions which are in the context specific models", "# reactions for reactions which are in the context specific models"];  
    title_plots = ["after discretization: ", "after mapping the gpr rules: ", "after applying " + project.models.(modelList(1)).settings.script_parameters.consensus_proportion + " consensus proportion"];  

    numDatasets = length(datasets);

    % Determine all unique discretization values across datasets (excluding 13)
    % the value 13 has been set to indicate that the rxn/gene is not in the
    % model, so the discretization is not shown, in these figures only the
    % discretization is shown of the genes/rxns in the model, the figures
    % for all genes, rxns are done in the QC script when preparing the data
    % for the model creation !!
    all_values = [];
    for k = 1:numDatasets
        all_values = union(all_values, setdiff(unique(datasets{k}), 13));
    end
    all_values = sort(all_values);  % e.g., [-1 0 1]
    
    % Define a colormap with one color per value
    cmap = lines(length(all_values));
    
    % Create figure
    plots.data_discretization = figure('Color','w','Position',[100 100 2000*numDatasets 2000],...
                                       'Visible','off');
    
    for k = 1:numDatasets
        dataset = datasets{k};
        xlabel_plot = xlabels_plots(k);
        xtick_plot = xticks_plots{k};
        ylabel_plot = ylabels_plots(k);
        title_plot = title_plots(k);
        
        % counts per category per sample (make sure all_values are included)
        counts = zeros(length(all_values), size(dataset,2));
        for i = 1:length(all_values)
            counts(i,:) = sum(dataset == all_values(i), 1);
        end
    
        % stacked barplot in subplot
        ax = subplot(1, numDatasets, k);
        b = bar(ax, counts', 'stacked');
    
        % Assign consistent colors
        for i = 1:length(all_values)
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
        xlabel(xlabel_plot,'FontSize',16)
        ylabel(ylabel_plot,'FontSize',16)
        title(title_plot,'FontSize',18)
    
        xticks(1:length(xtick_plot))
        xticklabels(regexprep(xtick_plot, "_", "-"))
    end
    
    % Single legend for the whole figure
    lgd = legend(string(all_values), 'Location','northeastoutside');
    lgd.FontSize = 14;
    lgd.Title.String = "Discretization status";


    % -- Visualization: Get the jaccard similarity on basis of the
    % gene,metabolite and reaction presence in the corresponding models
    % How similar are my models structuraly, which models are more similar
    % to each other than others ? 
   
    
    [rxn_presence,rxn_mapping] = getOrderedFeatureMatrix(project,modelList,"rxns", reference_model);
    structure_analysis.rxn_mapping_table = array2table(rxn_mapping,"VariableNames",modelList,"RowNames",string(project.models.(reference_model).model.rxns));
   

    for field_to_investigate = ["genes", "mets", "rxns"]
        [ordered_feature, ~] = getOrderedFeatureMatrix( ...
            project, modelList, field_to_investigate, reference_model);
    
        % Plot Venn / Heatmap of intersections based on presence
        plots.intersections.(field_to_investigate) =  plotFlexibleVenn( ...
                                                                    ordered_feature, ...
                                                                    structure_analysis.modelNames, ...
                                                                    "Structural model comparison: " + field_to_investigate + " presence");
    
        % get the jaccard distances - based on reaction presence
        % Compute Jaccard distances
        plots.jaccard_dist.(field_to_investigate) = figure('Position',[20 20 700 300],'Visible','off');
 
        Jacc_distance = 1 - squareform(pdist(ordered_feature','jaccard'));
        title_fig = "Jaccard similarity of " + field_to_investigate + " presence (0 or 1) between models";
        
        % Create heatmap
        h = heatmap(structure_analysis.modelNames, structure_analysis.modelNames, Jacc_distance);
        
        % Set font sizes
        h.FontSize = 20;           
        h.Title = title_fig;    
    end
    

    % -- Visualization: Rxns presence jaccard similarity per pathway ? 
   
    pathway_names = string(project.models.(reference_model).model.subSystems);
    
    pathway_wise_jaccard_sim = [];
    for x=unique(pathway_names)'
        rxn_presence_pathway = rxn_presence(find(matches(pathway_names,x)),:);
        Jacc_distance = 1 - pdist(rxn_presence_pathway','jaccard');
        pathway_wise_jaccard_sim = [pathway_wise_jaccard_sim,Jacc_distance];
    end
    
    
    % -- Visualization: Get reaction presence for each model in comparison
    % to the defined reference model -> visualization per subsystem
    % Where does the difference I see in the jaccard plot come from ? form
    % which subsystem, which subsystem is most different in pairwise
    % comparison ? 
        
    pathways = string(project.models.(reference_model).model.subSystems); % get pathways from reference model
    unique_pathways = unique(pathways); 

    % for every pathway get the matrix identifying the rnxs from reference
    % model in this pathway
    groups = arrayfun(@(x) find(pathways == x), unique_pathways, 'UniformOutput', false);
    num_groups = length(groups);
    G = zeros(num_groups, size(ordered_feature,1));

    for g = 1:num_groups
        G(g, groups{g}) = 1;
    end
    
    % get the presence per subsystem in the context specific models 
    % by mulitplying the rxns presence for each subsystem (matrix ordered features) 
    % with the matrices defining the subsystem for every rxns
    pathway_counts = array2table(G * ordered_feature, ...
                                 'VariableNames', structure_analysis.modelNames,...
                                 'RowNames',cellstr(unique_pathways));
    % add reference model count to be able to make a relative abundance
    pathway_counts.reference_model = groupcounts(pathways);
    
    relative_counts = array2table(pathway_counts{:,1:end-1} ./ pathway_counts.reference_model, ...
                                 'VariableNames', structure_analysis.modelNames,...
                                 'RowNames',cellstr(unique_pathways));
    
    % get the idx of the most variant pathways in terms of rxns presence
    relative_counts.row_var = var(relative_counts{:,:}, 0, 2);
    pathway_counts.row_var = var(pathway_counts{:,1:end-1},0,2);
    pathway_counts = pathway_counts(pathway_counts.reference_model < 1000,:);
    % Get indices of top n highest variance rows

    
    pathway_counts = sortrows(pathway_counts,"row_var","descend");
    pathway_counts = pathway_counts(find(pathway_counts.row_var ~= 0),:);
    relative_counts = relative_counts(pathway_counts.Properties.RowNames,:);
    % plot top 20 most variant pathways between the choosen models
    data = relative_counts{:,1:end-1};
    rowNames = string(relative_counts.Properties.RowNames);
    colNames = structure_analysis.modelNames;

    %%%%%%%%%%

    plots.reaction_pathway_presence = figure('Position',[20 20 700 300],'Visible','off');
    tiledlayout(1,4, ...
        'TileSpacing','compact', ...
        'Padding','compact')
    
    % ---- Bar plot (LEFT) ----
    ax1 = nexttile(1);
    barh(pathway_counts.reference_model)
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
    cmap = get_color_pallette();
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
            % (or multiply relative_counts by reference_model if needed)
            value = pathway_counts{i,j}; % +1 because first column is reference_model
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

    %%% Visualize the core reaction per model
    data = struct2cell(structfun(@(x) [ length(x.core_reactions) - sum(ismember(x.core_reactions, x.model.rxns)); ...
                                            sum(ismember(x.core_reactions, x.model.rxns));...
                                            length(x.model.rxns) - sum(ismember(x.core_reactions, x.model.rxns));...
                                            length(x.model.rxns)], ...
                                            models_list, 'UniformOutput', false));
    data = [data{:}];
    
    % ---- Create layout ----
    upper_data = data(2:3,:);

    %figure
    categories = fieldnames(models_list)';  % model names
    %t = tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact');

    % -- Visualization: get a sense of how many of the core reactions made
    % it into your model. In theory 100 percent of the reactions should be
    % in, but in practice this is not gona happen, adjusting the
    % thresholding/discretization in the preprocessing of the data for
    % rfastcormics could change the precentages you see in the following
    % figure
        

    
    plots.core_reactions = figure('Visible','off');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
    
    % --- first barplot
    ax1 = nexttile(1);
    bar(upper_data', 'stacked')
    
    
    % Labels
    set(gca, 'XTickLabel', categories, 'FontSize', 14)  % increase tick labels font
    xlabel('Models', 'FontSize', 14)
    ylabel('# rxns', 'FontSize', 14)
    
    % Legend
    legend({"non-core reactions", "core reactions"}, 'Location','northwest', 'FontSize', 14)
    
    % Title
    title('Core and non-core reactions per model', 'FontSize', 14)

    % ---- Compute percentages of upper stack ----
    total = sum(upper_data,1);                   % total per model
    percent_upper = 100 * upper_data(2,:) ./ total;  % percentage of upper bar
    
    % ---- Add text labels on top of upper bars ----
    for i = 1:size(upper_data,2)  % loop over models
        % x-position is the bar center, y-position is height of lower + upper
        x = i;
        y = upper_data(1,i) + upper_data(2,i)/2;  % middle of upper stack
        text(x, y, sprintf('%.1f%%', percent_upper(i)), ...
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
    percent_upper = 100 * data(2,:) ./ total;  % percentage of upper bar
    
    % ---- Add text labels on top of upper bars ----
    for i = 1:size(data,2)  % loop over models
        % x-position is the bar center, y-position is height of lower + upper
        x = i;
        y = data(1,i) + data(2,i)/2;  % middle of upper stack
        text(x, y, sprintf('%.1f%%', percent_upper(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold');
    end

    % Tile 3 (THIS IS THE KEY)
    ax3 = nexttile(3,[1 2]);   % column 2, span both rows
    axis(ax3,'off')
    hold(ax3,'on')
    

    core_reactions_included = struct2cell(structfun(@(x) x.core_reactions(find(ismember(x.core_reactions, x.model.rxns)))', ...
                                            models_list, 'UniformOutput', false));
    core_reactions_included = unique([core_reactions_included{:}]);
    
    core_presence = structure_analysis.rxn_mapping_table{core_reactions_included,:} ~= 0;
    [figV,idx_inter_outersections,~] = plotFlexibleVenn(core_presence, structure_analysis.modelNames, ... 
                                                        "Structural model comparison: core rxns presence");

    
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
        [figV,idx_inter_outersections,~] = plotFlexibleVenn( ...
            core_presence, structure_analysis.modelNames, ...
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

    
    if string(class(figV)) == 'matlab.ui.Figure'

        % -- Visualization: Looking in deeper into the core reactions, the core
        % is what is defined by the data, therefore portrays the underlying
        % biological chnages, so the question is which reactions are part of
        % the outer and intersections we saw in the previous venn/intersection
        % diagramm ? are the differences in core reactions only due to
        % exchange/import ? transporters ? This should be avoided!
            
        % create an upsetr plot for the all the inter and outersections
        % filter out the main intersection -> the one with the longest name
        names_intersections = fieldnames(idx_inter_outersections);
        [~,all_intersection] = max(cellfun(@(x) length(x), names_intersections));
        idx_inter_outersections = rmfield(idx_inter_outersections, names_intersections(all_intersection));
        % now get the pathway of every entry
        inter_outersections_pathways = structfun(@(x) string(project.models.(reference_model).model.subSystems(x)),...
                                                 idx_inter_outersections,'UniformOutput',false);
        C = struct2cell(inter_outersections_pathways);
        unique_pathways = unique(vertcat(C{:}));
        
    
        % Preprocess pathways: collapse transport
        S = structfun(@(x) regexprep(x,"^Transport.*","Transport"), inter_outersections_pathways, 'UniformOutput', false);
        pathways_unique = unique(regexprep(unique_pathways,"^Transport.*","Transport"));
        
        barNames = string(fieldnames(S));
        nBars = numel(barNames);
        
        % Build count matrix
        Y = cellfun(@(b) sum(S.(b)' == pathways_unique, 2)', barNames', 'UniformOutput', false);
        Y = cat(1, Y{:});
        
        % Sort bars by total counts (descending)
        [~, sortIdx] = sort(sum(Y,2), 'descend');
        Y = Y(sortIdx,:);
        barNames_sorted = barNames(sortIdx);
        
        % Plot
        plots.core_reactions_intersections = figure('Visible','off');
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
        ax.XTickLabel = regexprep(barNames_sorted, "_", " ");
        ax.FontSize = 20;
        xlabel('Model intersections/outersections','FontSize',20)
        ylabel('# Core Reactions','FontSize',20)
        title('Count of Core reactions per pathway and intersection/outersection','FontSize',20)
        
        lgd = legend(pathways_unique, 'Location','northeast');
        lgd.FontSize = 20;
    end

    % Further Visualizations ? #TODO ? 
    structure_analysis.plots = plots;
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
        % using rfastcormics function to map discretized data to the rxns
        mapping = mapExpressionToModel( ...
            modelStruct.model, ...
            modelStruct.discretized_data.values, ...
            modelStruct.settings.dico, ...
            string(modelStruct.discretized_data.gene_names), ...
            1);
        
        numberOfSamples = size(mapping, 2);
        % store it per sample, column in the discretized expression matrix
        modelStruct.mappedDiscRxns_sample = mapping;
        % and also as global mapping for the model, by multiplying it with
        % the consensus porportion used for the model construction 
        % parameters used in the model construction can be accessed in the 
        % settings. slot of each individual model
        
        % definition of initialCore reactions
        modelStruct.mappedDiscRxns = sum(mapping == 1, 2) >= (modelStruct.settings.script_parameters.consensus_proportion * numberOfSamples);
        % definition of the notExpressed genes
        notExpressed = find(sum(mapping == -1, 2) >= (modelStruct.settings.script_parameters.consensus_proportion * numberOfSamples));
        modelStruct.mappedDiscRxns = int32(modelStruct.mappedDiscRxns);
        modelStruct.mappedDiscRxns(notExpressed) = -1;
        % The definition of the unexpressed and initialCore rxns is done as
        % performed in rFASTCORMICS_v2
    end
end


function modelStruct = reorderDiscretizedToMatchGeneOrder(modelStruct)
    % This function reorders the discretized data slot of every model to
    % have the rows/genes in the same order as in the associated
    % model.genes slot. In addition it adds the gene symbol to the discretized data slot.
    % The discretized data is given back as a matrix. 
    if any(contains(fieldnames(modelStruct),"discretized_data"))
        if size(modelStruct.discretized_data,1) ~= length(string(modelStruct.settings.dico.SYMBOL))% in case this is the second time you run a modelComparison, the disc slot is already ordered so we skip this function in that case
            gene_symbol = string(modelStruct.settings.dico.SYMBOL);
            gene_id_in_model = string(modelStruct.settings.dico.gene_id_in_model);
            discTbl     = modelStruct.discretized_data;   % table with gene_names + data
            % Map discretized genes into full gene list
            [isPresent, idx] = ismember(gene_symbol, string(discTbl.gene_names));
            % Preallocate output table with NaNs
            outTbl = zeros(numel(gene_symbol), size(discTbl.values,2));
            % Fill rows that exist
            outTbl(isPresent,:) = discTbl{idx(isPresent), "values"};
            % Add gene names as first column
            modelStruct.discretized_data = table(gene_id_in_model,gene_symbol,outTbl, 'VariableNames',...
                                                 ["gene_id_in_model",string(modelStruct.discretized_data.Properties.VariableNames)]);
        end
    end
end



function fva_lower = getLowerTriangleBlock(fva_sim_rxns)
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
    % fva_sim_rxns: n x n cell array of comparisons
    % Returns a compact lower-triangle block cell array

    n = size(fva_sim_rxns,1);
    if n ~= size(fva_sim_rxns,2)
        error('Input must be a square cell array');
    end

    % Count how many rows/columns to keep: remove first row if needed
    % General: keep lower-triangle below the first row
    rowsToKeep = 2:n;       % always skip the first row
    colsToKeep = 1:(n-1);   % always skip the last column

    % Take the lower-triangle block
    fva_lower = fva_sim_rxns(rowsToKeep, colsToKeep);
end

function results = pathway_enrichment(sets, metric_matrix,feature_names)
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


    [metricSorted, sortIdx] = sort(metric_matrix,'descend');
    rxnsSorted = feature_names(sortIdx);
    N = numel(feature_names);
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


function Results = get_enrichment_table(project,modelList,fva_sim_rxns,idx_to_keep,reference_model, subSystems)
    % This function visualizes the enrichment results in a dotplot!!
    % #TODO: better documentation of the function1!!!
    arguments
        project
        modelList
        fva_sim_rxns 
        idx_to_keep
        reference_model
        subSystems =[]
    end
    if isempty(subSystems)
        subSystems = string(project.models.(reference_model).model.subSystems(idx_to_keep)); 
    end
    rxns       = string(project.models.(reference_model).model.rxns(idx_to_keep));        

    [uniqSubs, ~, idx] = unique(subSystems);

    subStruct = struct();

    for k = 1:numel(uniqSubs)
        subName = uniqSubs{k};
        fieldName = matlab.lang.makeValidName(subName);
    
        % Collect reactions belonging to this subsystem
        subStruct.(fieldName).name = subName;
        subStruct.(fieldName).rxns = rxns(idx == k);
    end
    
    fvaSim = getLowerTriangleBlock(fva_sim_rxns);
 
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
    
        if isempty(x) || isempty(y)
            continue
        end
        
        Results.(string(y)) = pathway_enrichment(subStruct , x(idx_to_keep),rxns);

    end
     
end

function dotplot(NES_tbl,FDR_tbl)
    % #TODO better documentation of the function!!


    % Extract pathways and comparisons
    pathways = regexprep(string(NES_tbl.Properties.RowNames), "_", " ");
    comparisons = string(NES_tbl.Properties.VariableNames);
    
    % Extract numeric matrices
    NES = NES_tbl{:,:};
    FDR = FDR_tbl{:,:};
    
    nP = numel(pathways);
    nC = numel(comparisons);

    % Create grid coordinates
    [X, Y] = meshgrid(1:nC, 1:nP);
    x = X(:);
    y = Y(:);
    
    nesVals = NES(:);
    fdrVals = FDR(:);
    
    % Dot size proportional to |NES|
    dotSize = abs(nesVals) * 200;   % increase for pronounced size differences
    
    % Prepare FDR for coloring
    cVals = fdrVals;            % base colors from FDR
    cVals(cVals > 0.05) = NaN;  % values >0.05 will be grey
    
    % Create figure
    figure('Position',[100 100 900 500])
    hold on
    
    % Scatter plot for FDR ≤ 0.05
    scatter(x(~isnan(cVals)), y(~isnan(cVals)), dotSize(~isnan(cVals)), cVals(~isnan(cVals)), 'filled')
    
    % Scatter grey dots for FDR > 0.05
    scatter(x(isnan(cVals)), y(isnan(cVals)), dotSize(isnan(cVals)), [0.7 0.7 0.7], 'filled')
    
    % Axes formatting
    xticks(1:nC)
    xticklabels(regexprep(comparisons,"_", " vs "))
    yticks(1:nP)
    yticklabels(pathways)
    
    xlabel('Model comparison')
    ylabel('Pathway')
    
    xlim([0.5, nC + 0.5])
    ylim([0.5, nP + 0.5])
    set(gca,'YDir','reverse','FontSize',18)
    title("Pathway enrichment (dot size = |NES|, color = FDR)")
    
    % Colorbar with red to blue
    nColors = 256;
    cmap = [linspace(1,0,nColors)' linspace(0,0,nColors)' linspace(0,1,nColors)']; 
    colormap(cmap)        % red (low) -> blue (high)
    caxis([0 0.05])       % fix color scaling
    cb = colorbar;
    cb.Label.String = 'FDR';
    cb.FontSize = 14;
    
    box on
end



function FVA_sim_values_hist(fva_sim_rxns, modelList)
    % This function visualizes the FVA values in a histogramm per
    % comparison. These histogramms give us an indication of how similary
    % models are in their FVA min and max boundaries, and although the
    % heatmap summing up the FVA values also hast the same intuition, this
    % visualization also enables to see what similarity values are mostly
    % occuring. Is the difference in the overall similarity mainly due to a
    % few reactions having very low values, or do we see a lot of mean fva
    % similarity values per reaction ? 
    % #TODO improve function documentation!!
    fva_lower2x2 = getLowerTriangleBlock(fva_sim_rxns);
    
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
    
    [nRows, nCols] = size(fva_lower2x2);
    
    % Create tiled layout
    t = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
    
    for i = 1:nRows
        for j = 1:nCols
            nexttile((i-1)*nCols + j)
    
            data = fva_lower2x2{i,j};

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
    sgt.FontSize = 20;
    

end

