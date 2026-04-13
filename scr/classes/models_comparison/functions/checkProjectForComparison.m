function project = checkProjectForComparison(project, modelList,referenceModel)

    % check that object contains models field 
    if ~ismember("models",fieldnames(project))
        error("Your struct is not an appropriate input for the fastcore pipeline cause it does not entail a .models slot!")
    end

    % check that object contains the specified models 
    if any(~ismember(modelList,fieldnames(project.models)))
        error("One or multiple of the models you defined in the modelList input are not listed in the project.models slot of the fastcore project struct. Check that you have given the correct model names as an input to modelList!")
    end
    
    % check that every model is a cobra layout model ? 
    for m = [modelList,referenceModel];
        errors = verifyModel(project.models.(m).model);
        if ismember("Errors", fieldnames(errors))
            error("The model defined in the slot of: " + m  + " does not have a proper model stored under project.models." + m + ".model! There has to be a proper cobra model stored in that slot!!")
        end
    end

    % check that the dico is the same length as the genes in the model 

    % check that the expression data has the same length in every model
    % slot -> needs to be in the same order 
    


end