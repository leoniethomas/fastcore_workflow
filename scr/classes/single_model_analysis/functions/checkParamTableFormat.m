function checkParamTableFormat(tbl)
% Checks that a table has columns Parameter, Analysis, Value
% and that these columns contain only text.
% Throws an error and stops if a problem is detected.

    % Check that input is a table
    if ~istable(tbl)
        error('Input must be a table.');
    end

    % Check number of columns
    expectedCols = {'Parameter', 'Analysis', 'Value'};

    % Check column names
    actualCols = tbl.Properties.VariableNames;
    for i = 1:length(expectedCols)
        if ~strcmp(actualCols{i}, expectedCols{i})
            error('Column %d: expected name ''%s'', found ''%s''.', ...
                  i, expectedCols{i}, actualCols{i});
        end
    end

    % Check that each column contains only text
    for i = 1:length(expectedCols)
        col = tbl.(expectedCols{i});
        if ~(iscellstr(col) || (isstring(col) && isvector(col)))
            error('Column ''%s'' must contain only text. Found type: %s.', ...
                  expectedCols{i}, class(col));
        end
    end

end