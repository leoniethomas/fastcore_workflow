function tbl = readInParamTable(filePath)
    % This function reads in the csv file entailing the defaultParameters. 
    % Use this in case you encounter a problem with the Format or content
    % of the ParameterTable. The Problem might be due to the file not
    % being read in properly. 

    opts = detectImportOptions(filePath);
    opts.VariableTypes = {'string','string', 'string'}; % making sure that the last column with the values is read in as a character 
    % Ensure all rows are read (no early stopping)
    opts.DataLines = [2 Inf];
    tbl = readtable(filePath, opts);
    colNames = tbl.Properties.VariableNames;
    tbl = varfun(@cellstr, tbl);
    tbl.Properties.VariableNames = colNames;
end
