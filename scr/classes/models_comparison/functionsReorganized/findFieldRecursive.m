function [values, paths] = findFieldRecursive(S, fieldName, currentPath)

    if nargin < 3
        currentPath = "";
    end

    values = [];
    paths = {};

    if ~isstruct(S)
        return
    end

    for i = 1:numel(S)

        fields = fieldnames(S(i));

        for j = 1:numel(fields)

            thisField = fields{j};
            value = S(i).(thisField);

            if currentPath == ""
                thisPath = string(thisField);
            else
                thisPath = currentPath + "." + thisField;
            end

            % Found requested field
            if strcmp(thisField, fieldName)
                values = [values; value];
                paths{end+1} = char(thisPath);
            end

            % Recurse into structs
            if isstruct(value)

                [v, p] = findFieldRecursive(value, fieldName, thisPath);

                values = [values; v];
                paths = [paths, p];

            % Recurse into cells containing structs
            elseif iscell(value)

                for k = 1:numel(value)

                    if isstruct(value{k})

                        cellPath = thisPath + "{" + k + "}";

                        [v, p] = findFieldRecursive( ...
                            value{k}, fieldName, cellPath);

                        values = [values; v];
                        paths = [paths, p];
                    end
                end
            end
        end
    end
end