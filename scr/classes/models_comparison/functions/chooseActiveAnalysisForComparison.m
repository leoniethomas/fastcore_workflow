function [project,analysisID] = chooseActiveAnalysisForComparison(project,modelList,loopless,analysisID)
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
    %   - modelList:    list of models for which a active analysis should
    %                   be defined
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
        
        analysis = mod.analysis.(analysisID(m));
        for slot=1:length(fieldnames(analysis))
            slot_names = string(fieldnames(analysis));
            project.models.(modelList(m)).analysis.(slot_names(slot)) = analysis.(slot_names(slot));
        end

        if loopless & ismember('sampling', fieldnames(analysis)) 
            project.models.(modelList(m)).analysis.sampling.samples = project.models.(modelList(m)).analysis.sampling.cycleFreeFlux.samples_ll;
        end
        
    end
end