function validateOptionalSettings(optSettings)
    accepted = {'func', 'medium', 'notMediumConstrained'};
    actual = fieldnames(optSettings);
    extra = setdiff(actual, accepted);
    if ~isempty(extra)
        error("Not accepted fields in optionalSettings: %s", strjoin(extra, ", "));
    end
end