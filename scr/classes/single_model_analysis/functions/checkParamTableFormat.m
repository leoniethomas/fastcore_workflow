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

    % Check that all the lines in the csv were read into the table 
    expectedNumberOfLinesInCSV = 30;
    if size(tbl,1) < expectedNumberOfLinesInCSV
        error('The read in Parameter Table needs to entail at least %s Rows. The current Table only entails %s Rows. Go back and make sure that when reading in the csv files you have as many rows in the read in table in your Matlab workspace as there are in the csv file. Check if using the function: readInParamTable.m when reading in your csv will fix the issue.', ...
                    string(expectedNumberOfLinesInCSV),string(size(tbl,1)));
    end

    % Check that each column contains only text
    for i = 1:length(expectedCols)
        col = tbl.(expectedCols{i});
        if ~(iscellstr(col) || (isstring(col) && isvector(col)))
            if expectedCols{i} == "Value" && tbl.Properties.VariableTypes{i} == "double"
                tbl.Value = string(tbl.Value);
            else
                error('Column ''%s'' must contain only text. Found type: %s.', ...
                    expectedCols{i}, class(col));
            end
        end
    end

end