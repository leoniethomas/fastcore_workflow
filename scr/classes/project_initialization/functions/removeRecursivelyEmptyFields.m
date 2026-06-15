function s = removeRecursivelyEmptyFields(s, exclude)
% Recursively removes empty fields from a struct.
% exclude (cell of strings): field names to skip entirely (not removed,
% not recursed into).

    if nargin < 2
        exclude = {};
    end

    fn = fieldnames(s);
    for i = numel(fn):-1:1
        if ismember(fn{i}, exclude)
            continue
        end
        val = s.(fn{i});
        if isstruct(val) && ~istable(val)
            s.(fn{i}) = removeRecursivelyEmptyFields(val, exclude);
            if isFieldEmpty(s.(fn{i}))
                s = rmfield(s, fn{i});
            end
        elseif isFieldEmpty(val)
            s = rmfield(s, fn{i});
        end
    end
end

function tf = isFieldEmpty(val)
    % Empty struct
    if isstruct(val)
        tf = isempty(fieldnames(val));
        return
    end
    % Empty string
    if isstring(val)
        tf = (strlength(val) == 0);
        return
    end
    % Default: empty array ([]), cell, table, etc.
    tf = isempty(val);
end