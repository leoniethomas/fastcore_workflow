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
    project.comparisons.(comparison_name) = modelStructureComparison(project,modelList,reference_model);
    project.comparisons.(comparison_name).reference_model = reference_model; 
    % run functional model comparison

end


function structure_analysis = modelStructureComparison(project, modelList,reference_model)
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
    models = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    models = structfun(@(x) x.model, models, 'UniformOutput',false);
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

    % get the presence of rxns per subsystem in relation to the reference model
    
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
    pathway_counts = pathway_counts(:,1:end-1);
    
    % get the idx of the most variant pathways in terms of rxns presence
    row_var = var(relative_counts{:,:}, 0, 2);
    % Get indices of top n highest variance rows
    [~, sortedIdx] = sort(row_var, 'descend');
    top_var_pathways = relative_counts(sortedIdx(1:20),:);
    
    % plot top 20 most variant pathways between the choosen models
    figure
    heatmap(structure_analysis.modelNames,...
            string(top_var_pathways.Properties.RowNames),...
            top_var_pathways{:,:})
    title("relative counts of subsystem rxn occurence/reference model");

    
                
end


function [ordered_feature_matrix,ordered_rxn_matrix_idx] = getOrderedFeatureMatrix(project,modelList,field_to_investigate,reference_model,replacement_value)
    % this function brings the features for multiple models (for example the rxns presence (0
    % or 1)) into the same order than the reference model specified. 
    % Input: 
    % - project:                project object which is the output of
    %                           single_model_analysis script
    % - modelList:              names from models from which to get the
    %                           feature precence
    % - field_to_investigate:   feature presence of interest as string
    %                           "genes", "mets", or "rxns"
    % - reference_model:        reference model that give the order after
    %                           which the choosen feature will be ordered 
    % - replacement_value:      which value to put into the matrix. Just a
    %                           1 indicating that the given feature is present 
    %                           or not (0) in each of the choosen models, 
    %                           or the actuall values of sampling,fba, or fva ?  
    % 
    % Output: 
    % - ordered_feature_matrix: Matrix storing the wanted features 
    %                           (presence of genes,mets,rxns,or actuall 
    %                           fva,fba,sampling values) in the order of 
    %                           the reference model.
    %                           dim nxm with n = features in reference model
    %                                   and  m = models to compare
    % 
    % - ordered_rxn_matrix_idx: Tabel with the actual positions the rxns are 
    %                           stored at in the individual models. 
    %                           dim nxm with n = rxns in reference model
    %                                   and  m = models to compare
    % 
    arguments
        project
        modelList
        field_to_investigate (1,1) string ="rxns"
        reference_model  = strings(0)
        replacement_value (1,1) =1
    end
    % by default we bring all the models in the same order as the original model
    if isempty(reference_model)
        reference_model = project.models.orig_model.model;
    else
        reference_model = project.models.(reference_model).model;
    end

    models = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    models = structfun(@(x) x.model, models, 'UniformOutput',false);

    ordered_feature_matrix = struct2array(structfun(@(x) getOrderedFeature(x,reference_model,field_to_investigate,replacement_value), ...
                                             models,'UniformOutput',false));
    
    ordered_rxn_matrix_idx = struct2array(structfun(@(x) getOrderedFeature(x,reference_model,"rxns","idx"), ...
                                             models,'UniformOutput',false));
end

function ordered_feature = getOrderedFeature(model,reference_model,field_to_investigate,replacement_value)
    % this gives you back a vector in the length of the defined reference
    % data field (for example the length of rxns in the reference model)
    % and then give a 1 if this feature (for example the rxns) is present
    % in the specified model and a 0 if it is not present in the specified
    % model.
    % this function brings the features for a single model(for example the rxns presence (0
    % or 1)) into the same order than the reference model specified. 
    % Input: 
    % - model:                  model from which to obtain the feature
    % - reference_model:        reference model that give the order after
    %                           which the choosen feature will be ordered 
    % - field_to_investigate:   feature presence of interest as string
    %                           "genes", "mets", or "rxns"
    % - replacement_value:      which value to put into the vector. Just a
    %                           1 indicating that the given feature is present 
    %                           or not (0) in each of the choosen models, 
    %                           or the actuall values of sampling,fba, or fva ?  
    % 
    % Output: 
    % - ordered_feature:        a vector storing the wanted features 
    %                           (presence of genes,mets,rxns,or actuall 
    %                           fva,fba,sampling values) 
    %                           in the order of the reference model.
    %                           length(ordered_feature) == n with 
    %                           n = length(model.(field_to_investigate))
    % 
    arguments
        model
        reference_model
        field_to_investigate (1,1) string ="rxns"
        replacement_value (1,1) =1
    end

    ordered_feature_idx = zeros(length(reference_model.(field_to_investigate)),1);
    [~,idx] = ismember(model.(field_to_investigate),reference_model.(field_to_investigate));
    ordered_feature_idx(idx) = 1:numel(model.(field_to_investigate));
    
    % now we have a matrix that stores the idx of the features of the
    % models on the position of where the feature is in the reference
    % model.
    ordered_feature = ordered_feature_idx;
    if isnumeric(replacement_value)
        ordered_feature(ordered_feature_idx >0) = replacement_value;
    elseif isstring(replacement_value) || ischar(replacement_value)
       %replacement_values = model.analysis.(replacement_value).v;
    end

end




