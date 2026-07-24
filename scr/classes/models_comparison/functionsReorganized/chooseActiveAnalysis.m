function [project, analysisIDs] = chooseActiveAnalysis(project, modelList, analysisIDs, overwriteActive)
    % This function needs to be run in preparation for the modelsComparison
    % function. For the loaded project object, multiple analysis with a
    % different set of parameters can be performed. Before going into the
    % comparison, an active analysis needs to be choosen out of all the
    % singleModelAnalysis on grounds of which the comparison is then made.
    % In order to do so this function moves the analysis performed one slot up
    % in the struct structure so it is directly in the analysis slot.
    % This way the downstream analysis can be performed without specifying
    % the exact analysis name for every single model in the comparison.
    %
    % Inputs:
    %   - project: the project (struct)
    %   - modelList: list of models for which an active analysis should
    %     be defined (cell or string array)
    %   - analysisIDs: the name of the analysis for each of the models to
    %     be set as the active analysis for the following comparison.
    %     When no analysisIDs is given the most recent analysis performed
    %     for each model will be set as active.
    %   - overwriteActive: when using {'all'} (default) then before adding the
    %     chosen analysis to the default slot, all the objects that are in
    %     there are deleted. In case you want to just replace (for example)
    %     the FBA in the current default then you can just overwrite the
    %     FBA slot without deleting (for example the sampling) by setting
    %     overwriteActive to {'FBA'}.
    %
    % Outputs:
    %   - project: a project with a defined active analysis
    %   - analysisIDs: get the analysisID used returned, in the same order
    %     as the modelList given

    arguments
        project
        modelList (1,:) cell
        analysisIDs (1,:) cell = {}
        overwriteActive (1,:) cell = {'all'}
    end

    % Convert inputs to string arrays for internal use
    modelList = string(modelList);
    overwriteActive = string(overwriteActive);
    if isempty(analysisIDs)
        analysisIDs = strings(1,0);
    else
        analysisIDs = string(analysisIDs);
    end

    % Validate overwriteActive
    if ismember("all", overwriteActive) && numel(overwriteActive) > 1
        error("overwriteActive cannot contain 'all' together with other field names.");
    end

    %% Check format of models/analysisIDs and defining active analyses
    if isempty(analysisIDs)
        checkProjectFormat(project, modelList)

        % Defining the most recent analysis as active in case no
        % analysisIDs were provided
        for m = 1:numel(modelList)
            modelShortcut = project.models.(modelList(m));

            if isfield(modelShortcut, 'analysis')
                analysisFields = string(fieldnames(modelShortcut.analysis));

                % Check that analysisFields is not empty
                if isempty(analysisFields)
                    error("There is no analysis entry for this model: '" + modelList(m) + "'. Run singleModelAnalysis function in order to create an analysis field for the specified models.");
                end

                % Search for slots created by singleModelAnalysis
                isA = startsWith(analysisFields, "analysis_");

                % Check that there is at least one analysis_ field
                if ~any(isA)
                    error("There is no analysis_ entry for this model: '" + modelList(m) + "'. Run singleModelAnalysis function in order to create an analysis field for the specified models.");
                end

                % Validate date format for analysis_ fields
                analysisFieldsA = analysisFields(isA);
                dateStrs = extractAfter(analysisFieldsA, "analysis_");
                isValidDate = false(size(dateStrs));

                for k = 1:numel(dateStrs)
                    try
                        datetime(dateStrs(k), 'InputFormat', 'yyyyMMdd_HHmm');
                        isValidDate(k) = true;
                    catch
                        % not a valid date — will be reported below
                    end
                end

                % Report invalid IDs
                if any(~isValidDate)
                    invalidIds = analysisFieldsA(~isValidDate);
                    fprintf("Warning: The following analysis IDs in model '%s' are not valid dates (yyyyMMdd_HHmm) and will be ignored: %s", ...
                        modelList(m), strjoin(invalidIds, ", "));
                end

                % Check that there is at least one valid date
                if ~any(isValidDate)
                    error("There is no valid analysis_ entry (yyyyMMdd_HHmm) for this model: '" + modelList(m) + "'. Run singleModelAnalysis function in order to create an analysis field for the specified models.");
                end

                % compute time difference between now and the time the
                % analysis were created (only for valid dates)
                timeDiff = NaN(size(analysisFields));
                validIdx = find(isA);
                validIdx = validIdx(isValidDate);
                timeDiff(validIdx) = minutes(datetime("now") - ...
                    datetime(dateStrs(isValidDate), 'InputFormat', 'yyyyMMdd_HHmm'));
                [~, idx] = min(timeDiff); % choose the most recently performed one

                analysisIDs(m) = analysisFields(idx);

            else
                error("There is no analysis field for this model: '" + modelList(m) + "'. Run singleModelAnalysis function in order to create an analysis field for the specified models.");
            end
        end

    else
        checkProjectFormat(project, modelList, analysisIDs)
    end

    %% Creating/updating active field
    for m = 1:numel(modelList)
        modelShortcut = project.models.(modelList(m));
        srcAnalysis = modelShortcut.analysis.(analysisIDs(m));
        srcFields = fieldnames(srcAnalysis);

        if ~isfield(modelShortcut.analysis, 'active') || isequal(overwriteActive, "all")
            % (Re)initializing active field — full replacement
            modelShortcut.analysis.active = struct();

            % Copy all fields from the selected analysis into active
            for k = 1:numel(srcFields)
                fieldName = srcFields{k};
                modelShortcut.analysis.active.(fieldName) = srcAnalysis.(fieldName);
                % Add analysisId to structs and tables only
                if isstruct(modelShortcut.analysis.active.(fieldName))
                    modelShortcut.analysis.active.(fieldName).analysisId = analysisIDs(m);
                end
            end

        else
            % Selective replacement — only overwrite specified fields
            for i = 1:numel(overwriteActive)
                fieldName = overwriteActive(i);

                % Check that the field exists in the source analysis
                if ~isfield(srcAnalysis, fieldName)
                    error("Field '%s' does not exist in analysis '%s' for model '%s'. Cannot overwrite.", ...
                        fieldName, analysisIDs(m), modelList(m));
                end

                % Report whether it's a new field or an overwrite
                if ~isfield(modelShortcut.analysis.active, fieldName)
                    warning("Field '%s' does not exist in active analysis for model '%s'. It will be added from analysis '%s'.", ...
                        fieldName, modelList(m), analysisIDs(m));
                else
                    fprintf("Overwriting '%s' in active analysis for model '%s' with analysis '%s'.", ...
                        fieldName, modelList(m), analysisIDs(m));
                end

                % Copy the field
                modelShortcut.analysis.active.(fieldName) = srcAnalysis.(fieldName);
                % Add analysisId to structs and tables only
                if isstruct(modelShortcut.analysis.active.(fieldName)) 
                    modelShortcut.analysis.active.(fieldName).analysisId = analysisIDs(m);
                end
            end

            % Merge parameters (if parameters was not explicitly overwritten)
            if ~ismember("parameters", overwriteActive)
                if isfield(srcAnalysis, 'parameters') && ...
                   isfield(modelShortcut.analysis.active, 'parameters')
                    % Determine which analyses are being replaced
                    parametersReplace = setdiff(overwriteActive, "parameters");
                    if ~isempty(parametersReplace)
                        oldParams = modelShortcut.analysis.active.parameters;
                        newParams = srcAnalysis.parameters;

                        % Keep rows from oldParams that are NOT being replaced
                        idxKeep = ~ismember(string(oldParams.Analysis), parametersReplace);
                        oldParams = oldParams(idxKeep, :);

                        % Take rows from newParams that ARE being replaced
                        idxAdd = ismember(string(newParams.Analysis), parametersReplace);
                        newParamsToAdd = newParams(idxAdd, :);

                        % Combine
                        modelShortcut.analysis.active.parameters = [oldParams; newParamsToAdd];
                    end
                end
            end
        end

        % Write back to project (MATLAB structs are pass-by-value)
        project.models.(modelList(m)) = modelShortcut;
    end

end