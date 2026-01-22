function project = modelsComparison(project, modelList,analyses)
    % This function runs a set of comparative analysis in order to see the
    % differences between 2 or more models. 
    % Inputs: Project, Names of the models, List of wanted analysis (all by
    % default)
    % Available analysis:
    %   
    % Output : project with comparison field
    arguments
        project
        modelList (1,:) string = strings(0)
        analyses  (1,:) string = ["modelStructureComparison"]
    end
    
    % if no list of genes is given just compare all the models in the
    % project.models slot
    if isempty(modelList)
        modelList = project.models;
    end

    comparison_name = join(modelList, "_vs_");
    project.comparisons.(comparison_name) = struct();

    % loop over all the functions/analyses you want to be applied to the
    % object 
    project.comparisons.(comparison_name).(analyses) = modelStructureComparison(project,modelList);


end


function project = modelStructureComparison(project, modelList)
    % The structure comparison is a function that compares the models
    % listed on structural differences. Structural differences in the
    % context of Fastcore can be defined as the set of reactions that are
    % kept when running fastcore. This means we check for the existence of
    % rxns, metabolites and genes in the model, and their overlap between
    % models. In order 
    % Inputs: Project
    %   
    % Output : project with analysis in the corresponding comparison field
    arguments
        project
        modelList
    end

    % extract models you want to compare from the project object
    models = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    models = structfun(@(x) x.model, models, 'UniformOutput',false)

    % get rid of all the unused genes in the model in order to 

    % silently process all condition models
    % convert to table
    data = struct2array(structfun(@(x) {numel(x.rxns); numel(x.mets); numel(x.genes)}, ...
                           models, 'UniformOutput', false))';
            
    T = array2table(data, ...
                    'VariableNames', {'count_reactions','count_metabolites','count_genes'}, ...
                    'RowNames', string(fieldnames(models))');
    disp(T)

    % get the intersection outersection table from all the models
    field_to_investigate = "genes"
    % get presence matrix - all the entries in the same order

    ordered_feature_matrix = getOrderedFeatureMatrix(project,modelList);
    
    if size(ordered_feature_matrix,2) < 5
        plotFlexibleVenn(ordered_feature_matrix, fieldnames(models), "Structural model comparison: " + field_to_investigate + " presence")
    end
    
end


function ordered_feature_matrix = getOrderedFeatureMatrix(project,modelList,field_to_investigate,replacement_value, reference_model)
    % this function brings the features (for example the rxns presence (0
    % or 1)) into the same order than the reference model specified. 
    % Output: A matrix nxm with n = # models to compare and m = # of
    % features in the reference model specified (per default the orig_model)
    arguments
        project
        modelList
        field_to_investigate (1,1) string ="rxns"
        replacement_value (1,1) =1
        reference_model  = strings(0)
    end
    % by default we bring all the models in the same order as the original model
    if isempty(reference_model)
        reference_model = project.models.orig_model.model;
    else
        reference_model = project.models.(reference_model).model;
    end

    models = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    models = structfun(@(x) x.model, models, 'UniformOutput',false);

    ordered_feature_matrix = struct2array(structfun(@(x) getOrderedFeature(x,reference_model,field_to_investigate), ...
                                             models,'UniformOutput',false));
end

function ordered_feature = getOrderedFeature(model,reference_model,field_to_investigate,replacement_value)
    % this gives you back a vector in the length of the defined reference
    % data field (for example the length of rxns in the reference model)
    % and then give a 1 if this feature (for example the rxns) is present
    % in the specified model and a 0 if it is not present in the specified
    % model.
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
    elseif isstring(replacement) || ischar(replacement)
        replacement_values = model.analysis.(replacement_value).v;

    end

end




