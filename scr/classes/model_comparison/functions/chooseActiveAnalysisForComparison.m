function [project,analysisID] = chooseActiveAnalysisForComparison(project, modelList, loopless, emptyDefault, analysisID)
    % This function needs to be run in preparation for the modelsComparison
    % function. For the loaded project object, multiple analysis with a
    % different set of parameters can be performed. Before going into the
    % comparison an active analysis needs to be choosen out of all the
    % singleModelAnalysis on grounds of which the comparison is then made.
    % In order to do so this function moves the analysis performed one slot up
    % in the struct structure so it is directly in the analysis slot. 
    % This way the downstream analysis can be performed without specifying
    % the exact analysis name for every single model in the comparison. 
    % Input: 
    %   - project:      the fastcore project
    %   - modelList:    list of models for which an active analysis should
    %                   be defined
    %   - loopless:     indicate whether the loopless sampling should be
    %                   used for the analysis or not (default = 1)
    %   - emptyDefault: when using 1 (default) then before adding the
    %                   choosen analysis to the default slot, all the
    %                   objects that are in there are deleted. In case you
    %                   want to just replace (for example) the FBA in the
    %                   current default then you can just overwrite the 
    %                   FBA slot without deleting (for example the sampling) 
    %                   by setting emptydefault to 0.
    %   - analysisID:   the name of the analysis for each of the models to
    %                   be set as the active analysis for the following
    %                   comparison. When no analysisID is given the most recent
    %                   analysis perfomed for each model will be set as active.
    % Output:
    %   - project:      a fastcore project with a defined active analysis
    %   - analysisID:   get the analysisID used returned, in the same order
    %                   as the modelList given
    %
    arguments
        project
        modelList (1,:) string
        loopless =1
        emptyDefault =1
        analysisID (1,:) string = []
    end


    % check the modelList - are the specified models exsistent ? 
    validModels = string(fieldnames(project.models));
    if ~all(ismember(modelList, validModels))
        invalid = modelList(~ismember(modelList, validModels));
        error("Invalid model name(s): %s. Valid models are: %s", ...
            strjoin(invalid,", "), ...
            strjoin(validModels,", "));
    end
    
    % check the choosen analysis to set as active
    if isempty(analysisID) 
        % when the argument is empty the most recently performed analysis
        % is choosen for each model
        for m=1:numel(modelList)
           analysis_fields = string(fieldnames(project.models.(modelList(m)).analysis));
            
            % search for slots created by singleModelAnalysis -> prefix
            % analysis_timestamp
            isA = startsWith(analysis_fields,"analysis_"); 
            timeDiff = NaN(size(analysis_fields));
            
            % compute the time difference between now and the time the
            % analysis were created
            timeDiff(isA) = minutes( ...
                datetime("now") - ...
                datetime(extractAfter(analysis_fields(isA),"analysis_"), ...
                "InputFormat","yyyyMMdd_HHmm") );
            [~,idx] = min(timeDiff); % choose the most recently performed one
            if all(isnan(timeDiff))
                error("There is no analysis object for this model: " + modelList(m) + ". Run singleModelAnalysis function in order to create an analysis object for the specified models!");
            end
            analysisID(m) = analysis_fields(idx);

        end
    elseif ~isempty(analysisID) && numel(analysisID) ~= numel(modelList)
        error("If analysisID is provided, it must have the same number of elements as modelList.")
    end
    
    % assign the specified analysis for the specified models to the default
    % slot 
    for m=1:numel(modelList)
        mod = project.models.(modelList(m));
        % check that analysisID is defined correctly 
        if ~ismember(analysisID(m),string(fieldnames(mod.analysis)))
            error("The specified analysis ID: " + analysisID(m) + " does not exsist in model: " + modelList(m))
        end
        % check if the default analysis slot is defined from previous runs
        % of chooseActiveAnalysis -> by default this will be removed
        % in case you want to keep some of the previous slots you need to
        % define emptyDefaul = 0
        if ismember("active",string(fieldnames(project.models.(modelList(m)).analysis))) && emptyDefault
            project.models.(modelList(m)).analysis = rmfield(project.models.(modelList(m)).analysis, "active");
        end
        
        analysis = mod.analysis.(analysisID(m));
        slot_names = string(fieldnames(analysis));
        for slot=1:length(slot_names)
            
            if slot_names(slot) == "parameters" && ~emptyDefault && ismember("active",string(fieldnames(project.models.(modelList(m)).analysis)))
                parametersReplace = setdiff(slot_names,slot_names(slot));
                oldParams = project.models.(modelList(m)).analysis.active.(slot_names(slot)); 
                newParams = analysis.(slot_names(slot));
                idxValueReplace = ismember(oldParams.Analysis,parametersReplace);
                if all(string(oldParams.Parameter(idxValueReplace)) == string(newParams.Parameter(idxValueReplace)))
                    oldParams.Value(idxValueReplace) = newParams.Value(idxValueReplace);
                    analysis.(slot_names(slot)) = oldParams;
                end
            end
            project.models.(modelList(m)).analysis.active.(slot_names(slot)) = analysis.(slot_names(slot));
            
        end

        if loopless
            if isfield(project.models.(modelList(m)).analysis.active,"sampling")
                if isfield(project.models.(modelList(m)).analysis.active.sampling,"cycleFreeFlux")
                    project.models.(modelList(m)).analysis.active.sampling.samples = project.models.(modelList(m)).analysis.active.sampling.cycleFreeFlux.samplesLl;
                else
                    project.models.(modelList(m)).analysis.active.sampling.samples = project.models.(modelList(m)).analysis.active.sampling.samples_loopless; % pb here no?
                end
            else
                error("There is no sampling slot defined in the active analysis slot! Check that there is a sampling analysis specified there!")
            end
        end
        
    end
end