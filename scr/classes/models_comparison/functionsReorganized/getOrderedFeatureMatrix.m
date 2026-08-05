function [orderedFeatureMatrix, orderedRxnMatrixIdx, columnLabels] = getOrderedFeatureMatrix(project, modelList, referenceModel, fieldToInvestigate, replacementValue)
    % this function brings the features for multiple models (for example the rxns presence (0
    % or 1)) into the same order than the reference model specified. 
    % Input: 
    % - project:                project object which is the output of
    %                           singleModelAnalysis script
    % - modelList:              names from models from which to get the
    %                           feature precence
    % - fieldToInvestigate:   feature presence of interest as string
    %                           "genes", "mets", or "rxns"
    % - referenceModel:        reference model that give the order after
    %                           which the choosen feature will be ordered 
    % - replacementValue:      which value to put into the matrix. Just a
    %                           1 indicating that the given feature is present 
    %                           or not (0) in each of the choosen models, 
    %                           or the actuall values of sampling, fba, or fva ?  
    % 
    % Output: 
    % - orderedFeatureMatrix: Matrix storing the wanted features 
    %                           (presence of genes,mets,rxns,or actuall 
    %                           fva,fba,sampling values) in the order of 
    %                           the reference model.
    %                           dim nxm with n = features in reference model
    %                                   and  m = models to compare
    % 
    % - orderedRxnMatrixIdx: Tabel with the actual positions the rxns are 
    %                           stored at in the individual models. 
    %                           dim nxm with n = rxns in reference model
    %                                   and  m = models to compare
    % - columnLabels:        shows to which of the models each solution belongs     
    % 
    arguments
        project
        modelList
        referenceModel
        fieldToInvestigate (1,1) string = "rxns"
        replacementValue (1,1) = 1
    end
    
    %% by default we bring all the models in the same order as the original model
    referenceModel = project.models.(referenceModel).model;

    models = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    %models = structfun(@(x) x.model, models, 'UniformOutput', false);

    orderedFeatureMatrix = structfun(@(x) getOrderedFeature(x, referenceModel, fieldToInvestigate, replacementValue), ...
                                             models, 'UniformOutput', false);

    sizeFeatureMatrix = structfun(@(x) size(x, 2), orderedFeatureMatrix);
    orderedFeatureMatrix = struct2array(orderedFeatureMatrix);
    columnLabels = repelem(modelList, sizeFeatureMatrix');
    
    orderedRxnMatrixIdx = struct2array(structfun(@(x) getOrderedFeature(x, referenceModel, fieldToInvestigate, "idx"), ...
                                             models, 'UniformOutput', false));
end

function orderedFeature = getOrderedFeature(model, referenceModel, fieldToInvestigate, replacementValue)
    % this gives you back a vector in the length of the defined reference
    % data field (for example the length of rxns in the reference model)
    % and then gives a 1 if this feature (for example the rxn) is present
    % in the specified model and a 0 if it is not present in the specified
    % model.
    % this function brings the features for a single model(for example the rxns presence (0
    % or 1)) into the same order than the reference model specified. 
    % Input: 
    % - model:                  model from which to obtain the feature
    % - reference_model:        reference model that give the order after
    %                           which the choosen feature will be ordered 
    % - fieldToInvestigate:   feature presence of interest as string
    %                           "genes", "mets", or "rxns"
    % - replacement_value:      which value to put into the vector. Just a
    %                           1 indicating that the given feature is present 
    %                           or not (0) in each of the choosen models, 
    %                           or the actuall values of sampling,fba, or fva ?  
    % 
    % Output: 
    % - orderedFeature:        a vector storing the wanted features 
    %                           (presence of genes,mets,rxns,or actuall 
    %                           fva,fba,sampling values) 
    %                           in the order of the reference model.
    %                           length(ordered_feature) == n with 
    %                           n = length(model.(field_to_investigate))
    % 
    arguments
        model
        referenceModel
        fieldToInvestigate (1,1) string = "rxns"
        replacementValue (1,1) = 1
    end
    
    if contains(string(replacementValue), "discretizedData") || contains(string(replacementValue), "mappedDiscretizedRxns")
        % do not replace with 0 when values are not there, cause 0 means something in the discretization
        orderedFeatureIdx = zeros(size(referenceModel.(fieldToInvestigate), 1), size(referenceModel.(fieldToInvestigate), 2)) + 13;
    else
        orderedFeatureIdx = zeros(size(referenceModel.(fieldToInvestigate), 1), size(referenceModel.(fieldToInvestigate), 2));
    end
    
    [~, idx] = ismember(model.model.(fieldToInvestigate), referenceModel.(fieldToInvestigate));
    orderedFeatureIdx(idx) = 1:numel(model.model.(fieldToInvestigate));
    
    % now we have a matrix that stores the idx of the features of the
    % models on the position of where the feature is in the reference
    % model.
    orderedFeature = orderedFeatureIdx;
    if isnumeric(replacementValue)
        orderedFeature(orderedFeatureIdx > 0) = replacementValue;
    elseif isstring(replacementValue) || ischar(replacementValue)
        if replacementValue ~= "idx"
            replacementValue = strsplit(replacementValue, ".");

            replacementValues = model.(replacementValue(1));
            
            for x = replacementValue(2:end)
                replacementValues = replacementValues.(x);
            end
            if size(replacementValues, 1) ~= length(model.model.(fieldToInvestigate))
                error("The size of the slot you have chosen (rxns, genes) does not agree with the replacement value slot. Are you sure the replacement slot you choose really contains information about the chosen slot (genes, rxns,mets) ?")
            end
            if size(replacementValues, 2) == 1
                orderedFeature(idx) = replacementValues;
            else
                orderedFeature = repmat(orderedFeature, 1, size(replacementValues, 2));
                orderedFeature(idx, :) = replacementValues;
            end
        end
    end
end