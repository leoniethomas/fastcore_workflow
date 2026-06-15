function validateManuallySetBoundaries(boundaries)
    accepted = {'closedImports', 'closedExports', 'unconstrainedImports', 'unconstrainedExports'};
    actual = fieldnames(boundaries);
    extra = setdiff(actual, accepted);
    if ~isempty(extra)
        error("Not accepted fields in manuallySetBoundaries: %s", strjoin(extra, ", "));
    end
end