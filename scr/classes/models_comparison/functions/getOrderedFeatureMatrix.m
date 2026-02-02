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
    %models = structfun(@(x) x.model, models, 'UniformOutput',false);

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
    if contains(string(replacement_value), "discretized") | contains(string(replacement_value), "mappedDiscRxns")
        % do not replace with 0 when values are not there, cause 0 means something in the discretization
        ordered_feature_idx = zeros(size(reference_model.(field_to_investigate),1),size(reference_model.(field_to_investigate),2)) + 13;
    else
        ordered_feature_idx = zeros(size(reference_model.(field_to_investigate),1),size(reference_model.(field_to_investigate),2));
    end
    [~,idx] = ismember(model.model.(field_to_investigate),reference_model.(field_to_investigate));
    ordered_feature_idx(idx) = 1:numel(model.model.(field_to_investigate));
    
    % now we have a matrix that stores the idx of the features of the
    % models on the position of where the feature is in the reference
    % model.
    ordered_feature = ordered_feature_idx;
    if isnumeric(replacement_value)
        ordered_feature(ordered_feature_idx >0) = replacement_value;
    elseif isstring(replacement_value) || ischar(replacement_value)
        if replacement_value ~= "idx"
            replacement_value = strsplit(replacement_value,".");
            replacement_values = model.(replacement_value(1));
            for x = replacement_value(2:end)
                replacement_values = replacement_values.(x);
            end
            if size(replacement_values,1) ~= length(model.model.(field_to_investigate))
                error("The size of the slot you have choosen (rxns, genes) does not aggree with the replacement value slot! Are you sure the replacement slot you choose really contains information about the choosen slot (genes, rxns,mets) ??")
            end
            if size(replacement_values,2) ==1
                ordered_feature(idx) = replacement_values;
            else
                ordered_feature = repmat(ordered_feature, 1, size(replacement_values,2));
                ordered_feature(idx,:) = replacement_values;
            end
        end
    end

end