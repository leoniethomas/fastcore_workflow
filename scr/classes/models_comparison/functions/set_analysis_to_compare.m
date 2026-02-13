function project = set_analysis_to_compare(project,model_name,analysis_id)
    % This function sets a specified set of analysis to be the default
    % analysis, to be used in the downstream comparison pipeline. 
    % By default it uses the most recent analysis found in the
    % model.analysis slot 
    arguments
        project
        model_name (1,:) string
        analysis_id (1,:) string = []
    end
    if isempty(analysis_id)
        % search for youngest analysis in each model.analysis slot specified
        for model=1:numel(model_name)
           analysis_fields = fieldnames(project.models.(model_name(model)).analysis);
           analysis_fields = string(analysis_fields);

            isA = startsWith(analysis_fields,"analysis_");
            timeDiff = NaN(size(analysis_fields));
            
            timeDiff(isA) = minutes( ...
                datetime("now") - ...
                datetime(extractAfter(analysis_fields(isA),"analysis_"), ...
                "InputFormat","yyyyMMdd_HHmm") );
            [~,idx] = min(timeDiff);
            analysis_id(model) = analysis_fields(idx);

        end
    end

    for model=1:numel(model_name)
        
        analysis = project.models.(model_name(model)).analysis.(analysis_id(model));
        for slot=1:length(analysis_id)
            slot_names = string(fieldnames(analysis));
            project.models.(model_name(model)).analysis.(slot_names(slot)) = analysis.(slot_names(slot));
        end
    end
end