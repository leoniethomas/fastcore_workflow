function project = modelsComparison(project, modelList,reference_model,analyses,identifier)
    % This function runs a set of analysis for the comparison of the
    % specified genes.
    % A number of analysis are run: 
    % - structural analysis: based on the differential presence of
    %   metabolites, genes and reactions in the different models 
    % - functional analysis: based on the quantitative values like FVA,
    %   FBA, sampling 
    % 
    % Inputs: 
    %   - project:          the object which is the output of the single_model_analysis
    %                       entailing the results of fba,fva,sampling, single gene
    %                       deletion etc. for a single model 
    %   - modelList:        the list of Model names to be included in the comparison 
    %   - reference_model:  the reference model used to compute the relative feature presence to 
    %   - analyses:         the list of analyses which should be performed 
    %   - identifier:       a string, will be added as a postfix to the analysis name, can be choosen freely
    %                       default)
    %   
    % Output : 
    %   - project:          project object with a added comparison field entailing
    %                       all the output, modelcomparison information

    arguments
        project
        modelList (1,:) string = strings(0)
        reference_model (1,:) string = "orig_model"
        analyses  (1,:) string = ["modelStructureComparison"]
        identifier (1,:) string = string(datetime('now','Format','_yyyyMMdd_HHmmss'))
    end
    
    % if no list of genes is given just compare all the models in the
    % project.models slot
    if isempty(modelList)
        modelList = project.models;
    end
    
    % give the comparison the name of all models + a identifier choosen
    comparison_name = join(modelList, "_vs_") + identifier;

    % run structural model comparison
    project.comparisons.(comparison_name) = modelStructuralComparison(project,modelList,reference_model);
    project.comparisons.(comparison_name).reference_model = reference_model; 
    % run functional model comparison


    modelFunctionalComparison(project, comparison_name,analyses);
    

end

function modelFunctionalComparison(project, comparison_name,analyses)
    % The structure comparison is a function that compares the models
    % listed on structural differences. Structural differences in the
    % context of Fastcore can be defined as the set of reactions that are
    % kept when running fastcore. This means we check for the existence of
    % rxns, metabolites and genes in the model, and their overlap between
    % models.

    arguments
        project
        comparison_name
        analyses
    end

    modelList = project.comparisons.(comparison_name).modelNames;
    reference_model = project.comparisons.(comparison_name).reference_model;
    % first visualize the fba solution 
    % what are the variables I need from the fba solution 
    
    fba_objective_values = cell2mat(cellfun(@(x) project.models.(x).analysis.FBA.f(1,1) ,modelList,"UniformOutput",false));
    get_exchange_rxns_idx = find(findExcRxns(project.models.(reference_model).model));    
    
    figure
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

    %%% barplot for the biomass/defined objective 

    nexttile(1);
    bar(fba_objective_values)
    title('Model comparison: flux of optimized reaction')
    ylabel('Reaction flux value for objective function [mMol/(gDW*h)]')
    xlabel('Model')
    xticklabels(modelList)
    
    % import rxns -> upper limit applied to function
    get_flux_plot(project, comparison_name,get_exchange_rxns_idx, ...
                  'threshold_flux','upper','FVA',true,'reducedCost',true)
    % expor rxns -> upper limit applied to function
    get_flux_plot(project, comparison_name,get_exchange_rxns_idx,...
                  'threshold_flux','lower','FVA',true,'reducedCost',true)
    

    % compute Fluxvariability display in heatmap 
    %%% -> model wide 

    %%% -> per subsystem

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
    % models.

    arguments
        project
        modelList
        reference_model
    end

    % extract models you want to compare from the project object
    models_list = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    models = structfun(@(x) x.model, models_list, 'UniformOutput',false);
    structure_analysis.modelNames = string(fieldnames(models));


    % get the size of the different models
    data = struct2array(structfun(@(x) {numel(x.rxns); numel(x.mets); numel(x.genes)}, ...
                           models, 'UniformOutput', false))';
    array2table(data, ...
                    'VariableNames', {'count_reactions','count_metabolites','count_genes'}, ...
                    'RowNames', string(fieldnames(models))')
   
    %%%% get the indices of all the rxns from the different models
    %%%% specified to put them into one table

    
    [~,rxn_mapping] = getOrderedFeatureMatrix(project,modelList,"rxns", reference_model);
    structure_analysis.rxn_mapping_table = array2table(rxn_mapping,"VariableNames",modelList,"RowNames",string(project.models.(reference_model).model.rxns))
    
    %%%%%%%%%%%%%%%%%%%% get the intersections/outersection 
    %%%%%% similarities between models based on rxns,genes, mets

    for field_to_investigate = ["genes", "mets", "rxns"]
        [ordered_feature, ~] = getOrderedFeatureMatrix( ...
            project, modelList, field_to_investigate, reference_model);
    
        % Plot Venn / Heatmap of intersections based on presence
        plotFlexibleVenn( ...
            ordered_feature, ...
            structure_analysis.modelNames, ...
            "Structural model comparison: " + field_to_investigate + " presence");
        % plot Venn/Heatmap of the overall intersection/outersection of the
        % models based on rxn prescence 
        plotFlexibleVenn(ordered_feature,...
                             structure_analysis.modelNames, ... 
                             "Structural model comparison: " + field_to_investigate + " presence");
    
        % get the jaccard distances - based on reaction presence
        Jacc_distance = 1 - squareform(pdist(ordered_feature','jaccard'));
        title_fig = "Jaccard similarity of " + field_to_investigate + " presence (0 or 1) between models ";
        figure
        heatmap(structure_analysis.modelNames,structure_analysis.modelNames, Jacc_distance);
        title(title_fig);
    end
    
    %%%%%%%%%%%%%%%%%%%%% get presence of rxns in each subsytem/pathway
        
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

    figure
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
    ax2 = nexttile(2,[1 3]);
    
    imagesc(data)
    nColors = 256;
    whiteToBlue = [linspace(1,0,nColors)', linspace(1,0,nColors)', ones(nColors,1)];
    colormap(ax2, whiteToBlue)
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

    figure
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
    figV = plotFlexibleVenn(core_presence, structure_analysis.modelNames, ... 
                     "Structural model comparison: core rxns presence");

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


end







function essential_pathway_metabolites = get_essential_pathway_metabolites(pathway,project,reference_model)
    % source for pathway metabolite composition: https://github.com/sysbiolux/MidbrainOrganoid_Miro1_scMetabMod/blob/main/ModelAnalysis/compareFBA_v3_PAPER.m
    
    arguments
        pathway (1,1) string {mustBeMember(pathway,["Glycolysis","PPP","TCA_cytoplasm_only","Tricarboxylic_acid_cycle","Oxidative_phosphorylation","Fatty_acid_oxidation","Cholesterol","Dopamine_Tyrosine_metabolism","Fatty_acid_oxidation__all__crn","Fatty_acid_oxidation__all__coa"])}
        project
        reference_model
    end
        pathways.Glycolysis = ["glc_D[c]","g6p[c]","f6p[c]","fdp[c]","dhap[c]","g3p[c]","13dpg[c]","3pg[c]","2pg[c]","pep[c]","pyr[c]"];
        pathways.PPP=["ru5p_D[c]","s7p[c]","e4p[c]"];
        pathways.TCA_cytoplasm_only= ["cit[c]","icit[c]","akg[c]","succ[c]","fum[c]","oaa[c]"];
        pathways.Tricarboxylic_acid_cycle=["cit[m]","icit[m]","akg[m]","succoa[m]","succ[m]","fum[m]","mal_L[m]","oaa[m]"];
        pathways.Oxidative_phosphorylation=["nadh[m]","fadh2[m]","focytC[m]","q10h2[m]","atp[m]"];
        pathways.Fatty_acid_oxidation=["accoa[c]","accoa[m]","coa[c]","coa[m]"];
        %%
        temp=find(ismember(string(project.models.(reference_model).model.subSystems),'Fatty acid oxidation'));
        metList=findMetsFromRxns(project.models.(reference_model).model, ...
                                 project.models.(reference_model).model.rxns(temp))';
        temp=find(contains(metList,'coa'));
        pathways.Fatty_acid_oxidation__all__coa=string(metList(temp));
        %%
        temp=find(ismember(string(project.models.(reference_model).model.subSystems),'Fatty acid oxidation'));
        metList=findMetsFromRxns(project.models.(reference_model).model, ...
                                 project.models.(reference_model).model.rxns(temp))';
        temp=find(contains(metList,'crn'));
        pathways.Fatty_acid_oxidation__all__crn=string(metList(temp));
        %%
        pathways.Cholesterol=["sql[r]","zymst[r]","chsterol[r]","chsterol[c]"];
        pathways.Dopamine_Tyrosine_metabolism=["phe_L[c]","tyr_L[c]","34dhphe[c]","dopa[c]","34dhpac[c]","34dhpha[c]","34dhpe[c]","tym[c]"];

        essential_pathway_metabolites = pathways.(pathway);
end




