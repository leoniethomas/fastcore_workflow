function project = singleModelAnalysis(project, parameterTable, modelList, analyses, saveCheckpoint, resumeFromCheckpoint)
% This function runs the analysis on a model or a list or models and stores
% the results in a structure.
% Inputs: Project, Parameter Table (at least the default one), Names of the models, List of wanted analysis (all by
% default)
% Available analysis:
% FBA
% FVA
% Shadow prices % needs to be added
% Sampling
% Loopless sampling
% FDR correction for sampling
% SingleGeneDeletion
% DoubleGeneDeletion
% Enrichment % needs to be added
% Output : project with an analysis field
%% Arguments block
arguments
    project struct
    parameterTable table
    modelList {mustBeText} = {}
    analyses cell {mustBeVectorOrEmpty} = {}
    saveCheckpoint logical = true
    resumeFromCheckpoint logical = false
end

%% Check the structure of project
modelList = checkStructureForSingleModelAnalysis(project, modelList); % get valid models

%% Check format of parameter table
checkParamTableFormat(parameterTable)

%% Check that given analysis are part of the list
validAnalyses = {'FBA', 'FVA', 'sampling', 'loopless', 'kld', 'singleGeneDeletion', 'doubleGeneDeletion'};

if isempty(analyses) || (ischar(analyses) && isempty(analyses))
    toPerform = validAnalyses;
else
    isImplemented = ismember(analyses, validAnalyses);
    if ~all(isImplemented)
        notImplemented = analyses(~isImplemented);
        fprintf('The following analyses are not implemented and will not be performed:');
        for i = 1:length(notImplemented)
            fprintf(' - %s', notImplemented{i});
        end
    end
    toPerform = analyses(isImplemented);
end

% check that the loopless ID given are valid analysisIDs 
% check that the given analysis IDs for the loopless are actually
% available in the corresponding models 
if any("loopless" == analyses) && ~any("sampling" == analyses)
    % check that we have as many samplingIDs given as models
    analysisIDSampling = parameterTable.Value{parameterTable.Parameter == "samplingToUse" & parameterTable.Analysis == "loopless"};
    analysisID = split(string(analysisIDSampling),"/");
    assert(length(analysisID) == length(modelList));
    % check that the sampling analysisIDs are in the same order as the
    % models
    for m = 1:numel(modelList)
        mod = string(modelList(m));
        anaID = analysisID(m);
        if isfield(project.models.(mod),"analysis")
            if ~ismember(anaID, fieldnames(project.models.(mod).analysis))
                error("The analysisID (%s) for model %s can not be found! Check again the analysis IDs (samplingToUse parameter) you specified in the defaultParametersAnalysis.csv. The analysisIDs should be in the same order as their corresponding models are given in the modelList.",...
                      anaID,mod)
            end
        end
    end
end

% Checkpoint: optionally resume from last saved model
checkpointFile = 'singleModelAnalysis_checkpoint.mat';
startIdx = 1;

if resumeFromCheckpoint && exist(checkpointFile, 'file')
    loaded = load(checkpointFile);
    project = loaded.project;
    startIdx = loaded.i + 1;
    if startIdx <= numel(modelList)
        fprintf('Checkpoint loaded — resuming at model %d/%d (%s)\n', ...
                startIdx, numel(modelList), modelList{startIdx});
    end
end

%% Run analysis
for i = startIdx:numel(modelList)
    name = modelList{i};
    fprintf('Analysis running for model %s (%d/%d)\n', name, i, numel(modelList));

    if ~isfield(project.models.(name), 'analysis')
        project.models.(name).analysis = struct();
    end

    % Analysis id
    id = ['analysis_' char(datetime("now", "Format", "yyyyMMdd_HHmm"))];
    project.models.(name).analysis.(id) = struct();
    
    % select analysisID used for sampling in case loopless is run 
    if any("loopless" == analyses) & ~any("sampling" == analyses)
        parameterTable.Value{parameterTable.Parameter == "samplingToUse" & parameterTable.Analysis == "loopless"} = analysisID(i);
    else % if the sampling is run in the same analysis run, that sampling will be automatically used for the loopless sampling
        parameterTable.Value{parameterTable.Parameter == "samplingToUse" & parameterTable.Analysis == "loopless"} = id;
    end
    % Store the settings
    project.models.(name).analysis.(id).parameters = parameterTable;

    % Checkpoint: try/catch to save state before crash propagates
    try
        project = performAnalysis(project, parameterTable, name, toPerform, id);
    catch ME
        fprintf('Error on model %s : %s', name, ME.message);
        fprintf('Stack:');
        for s = 1:numel(ME.stack)
            fprintf('  %s (line %d)', ME.stack(s).name, ME.stack(s).line);
        end
        fprintf('Last valid checkpoint: model %d (%s).\n', i-1, modelList{i-1});
        fprintf('Partial project returned. Relaunch with resumeFromCheckpoint=true to resume.\n');
        return;
    end

    % Checkpoint: save after each successful model
    if saveCheckpoint
        save(checkpointFile, 'project', 'i', '-v7.3');
        fprintf('Checkpoint saved (model %d/%d)\n', i, numel(modelList));
    end
    % End checkpoint save

end

%% Checkpoint: clean up if all models completed
if saveCheckpoint && exist(checkpointFile, 'file')
    delete(checkpointFile);
    fprintf('All models analyzed — checkpoint deleted.');
end

end