function validateOptionalSettings(optSettings)
    accepted = {'func', 'medium', 'notMediumConstrained'};
    actual = fieldnames(optSettings);
    extra = setdiff(actual, accepted);
    if ~isempty(extra)
        error("Not accepted fields in optionalSettings: %s", strjoin(extra, ", "));
    end
    for i = 1:numel(accepted)
        if isfield(optSettings, accepted{i}) && ~iscell(optSettings.(accepted{i}))
            error('optSettings.%s must be a cell array.', accepted{i});
        end
    end
end