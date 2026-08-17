function activeAnalysisTable = getActiveAnalysisIDTable(project,modelList)

    % filter models once
    filteredModels = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    modelNames     = string(fieldnames(filteredModels));
    
    % get analysis fields excluding parameters
    analysisFields = setdiff(fieldnames(structfun(@(x) x.analysis.active, filteredModels, 'UniformOutput', true)), 'parameters');
    
    % build table
    activeAnalysisTable = struct();
    for f = analysisFields'
        activeAnalysisTable.(string(f)) = structfun(@(x) x.analysis.active.(string(f)).analysisId, filteredModels, 'UniformOutput', true);
    end
    
    % convert and transpose
    activeAnalysisTable = rows2vars(struct2table(activeAnalysisTable, 'RowNames', modelNames));
    activeAnalysisTable.Properties.RowNames = activeAnalysisTable.OriginalVariableNames;
    activeAnalysisTable = removevars(activeAnalysisTable, 'OriginalVariableNames');

end