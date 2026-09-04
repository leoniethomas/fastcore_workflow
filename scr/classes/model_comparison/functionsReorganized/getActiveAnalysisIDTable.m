function activeAnalysisTable = getActiveAnalysisIDTable(project,modelList)

    % filter models once
    filteredModels = rmfield(project.models, setdiff(fieldnames(project.models), modelList));
    modelNames     = string(fieldnames(filteredModels));
    
    % get analysis fields excluding parameters
    analysisInModels = structfun(@(x) fieldnames(x.analysis.active), filteredModels, 'UniformOutput', false);
    C = struct2cell(analysisInModels);
    analysisFields = C{1};
    for k = 2:numel(C)
        analysisFields = intersect(analysisFields, C{k});
    end
    analysisFields = setdiff(analysisFields,'parameters');
    
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