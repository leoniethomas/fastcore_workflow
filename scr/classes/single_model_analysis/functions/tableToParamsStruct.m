function params = tableToParamsStruct(defaultParametersAnalysis, analysisName, model)

    % Filter rows for the requested analysis
    rows = defaultParametersAnalysis( ...
        strcmp(string(defaultParametersAnalysis.Analysis), string(analysisName)), :);

    if isempty(rows)
        error('No parameters found for analysis "%s".', analysisName);
    end

    % Initialize output structure
    params = struct();

    % Loop over parameters
    for i = 1:height(rows)

        p = strtrim(rows.Parameter{i});   % Parameter name (possibly with dots)
        val = strtrim(rows.Value{i});    % Parameter value (string)

        % --- Parse value ---
        value = parseValue(val, model);

        % --- Handle nested fields like "options.nPointsReturned" ---
        fields = strsplit(p, '.');

        % Assign dynamically into nested struct
        params = assignNested(params, fields, value);
    end

    % -------- Nested helper functions --------

    function out = parseValue(v, model)
        % Try numeric
        num = str2double(v);
        if ~isnan(num)
            out = num;
            return
        end

        % Try to evaluate expressions involving model (e.g. size(model.S,2))
        if nargin > 1 && contains(v, "model")
            try
                out = eval(v);
                return
            catch
                % Fall through to string if eval fails
            end
        end

        % Otherwise keep as string
        out = v;
    end

    function s = assignNested(s, fields, value)
        % Recursive nested assignment
        if numel(fields) == 1
            s.(fields{1}) = value;
        else
            f = fields{1};

            if ~isfield(s, f) || ~isstruct(s.(f))
                s.(f) = struct();
            end

            s.(f) = assignNested(s.(f), fields(2:end), value);
        end
    end

end
