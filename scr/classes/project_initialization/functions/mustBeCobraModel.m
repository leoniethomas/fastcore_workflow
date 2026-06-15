function mustBeCobraModel(model)
    required = {'S','rxns','mets','lb','ub'};
    missing = required(~isfield(model, required));
    if ~isempty(missing)
        error("Your model is not a valid Cobra model. Missing fields: " + strjoin(missing, ", "));
    end
end