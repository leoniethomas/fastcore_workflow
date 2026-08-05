function [project, comparisonName] = modelsComparison(project, modelList, referenceModel, identifier, analyses)
% MODELSCOMPARISON Runs a set of analyses for the comparison of the specified models.
%
% A number of analyses are run:
% - structural analysis: based on the differential presence of
%   metabolites, genes and reactions in the different models
% - functional analysis: based on the quantitative values like FVA, FBA
% - sampling analysis: comparative analysis based on the sampling results
%
% Inputs:
% - project: the object which is the output of the single_model_analysis
%   entailing the results of fba, fva, sampling, single gene
%   deletion etc. for a single model
% - modelList: the list of Model names to be included in the comparison
% - referenceModel: the reference model used to compute the relative reaction presence
% - analyses: the list of analyses which should be performed
%   + structuralComparison: investigates the differences between the models
%     on a structural level, is a gene, metabolite, or rxn present or not?
%   + functionalComparison: investigates the differences between models on a
%     functional level, how much flux do reactions carry in FVA, FBA solutions?
%   + samplingComparison: investigates the differences between models sampling
%     solutions, investigates samples solution space
% - identifier: a string, will be added as a postfix to the analysis name,
%   can be chosen freely (default: current timestamp)
%
% Output:
% - project: project object with an added comparison field entailing
%   all the output, modelcomparison information
% - comparisonName: gives back the name of the comparison added

arguments
    project        (1,1) struct
    modelList      (1,:) string
    referenceModel (1,1) string
    identifier     (1,1) string = string(datetime('now', 'Format', '_yyyyMMdd_HHmmss'))
    analyses       (1,:) string {mustBeMember(analyses, {'structuralComparison', 'functionalComparison', 'samplingComparison', 'IDAREoutput'})} = "structuralComparison"
end

%% Check project and models format for comparison

% Check first the reference model
checkProjectFormat(project, referenceModel)

% Check then models and analysis to compare
% active analysis of each model will be used for functionalComparison
checkProjectFormat(project, modelList, repmat({'active'}, 1, numel(modelList)))

%% Reorder models according to their order of appearance in project.models
order = string(fieldnames(project.models));
[~, idx] = ismember(modelList, order);
[~, sortIdx] = sort(idx);

modelListOrdered = modelList(sortIdx);

%% Give the comparison the name of all compared models associated with a given identifier
comparisonName = join(modelListOrdered, "_vs_") + "__" + identifier;

%% Create comparisons slot if not already existing
if ~isfield(project, 'comparisons')
    project.comparisons = struct();
elseif isfield(project, 'comparisons') && ~isstruct(project.comparisons)
    warning('project.comparisons exists but is not a struct (current class: %s).', ...
            class(project.comparisons));
    reply = '';
    while ~ismember(lower(strtrim(reply)), {'y', 'n', 'yes', 'no'})
        reply = input('Delete project.comparisons and reinitialize as an empty struct? [y/n] ', 's');
    end
    if startsWith(lower(strtrim(reply)), 'y')
        project.comparisons = struct();
        disp('project.comparisons reinitialized as an empty struct.');
    else
        error('Operation cancelled by user. project.comparisons left unchanged (class: %s).', ...
              class(project.comparisons));
    end
end

%% Comparison
% Structural comparison is needed for functional and sampling comparisons,
% therefore we check whether it was already run, to not run
% it again if it was already created

if isfield(project.comparisons, comparisonName)

    % TO ADD: checking of the format
    % If format correct, it means there is a field called
    % "referenceModel" which is part of project.comparisons

    if isequal(referenceModel, project.comparisons.(comparisonName).referenceModel)
        if isfield(project.comparisons.(comparisonName), 'structuralAnalysisStatus') && (project.comparisons.(comparisonName).structuralAnalysisStatus == 1)
            % Structural analysis already run
            disp("Structural comparison already run.");
            % Perform the other analyses if specified in analyses

            % Functional comparison
            if any(matches(analyses, "functionalComparison"))
                disp("Running functional comparison.");
                project.comparisons.(comparisonName).functionalComparison.plots = functionalComparison(project, comparisonName);
            end

            % Sampling comparison
            if any(matches(analyses, "samplingComparison"))
                disp("Running sampling comparison.");
                project = samplingComparison(project, comparisonName);
            end
        else
            % Structural analysis not run yet, mandatory
            project = runAllComparisons(project, modelListOrdered, referenceModel, comparisonName, analyses);
        end
    else
        % WARNING: field exists with a different referenceModel
        warning('Comparison "%s" already exists with referenceModel = "%s" (requested: "%s").', ...
                comparisonName, project.comparisons.(comparisonName).referenceModel, referenceModel);

        reply = '';
        while ~ismember(lower(strtrim(reply)), {'y', 'n', 'yes', 'no'})
            reply = input('Overwrite the existing comparison with the new reference model? [y/n] ', 's');
        end

        if startsWith(lower(strtrim(reply)), 'y')
            % YES: overwrite and rerun everything
            project = runAllComparisons(project, modelListOrdered, referenceModel, comparisonName, analyses);
        else
            % NO: stop, user must choose a new identifier
            error('Operation cancelled by user. Comparison "%s" keeps referenceModel = "%s".Choose a different identifier to create a new comparison.', ...
                  comparisonName, project.comparisons.(comparisonName).referenceModel);
        end
    end

else
    % Comparison does not exist yet: initialize and run everything
    project = runAllComparisons(project, modelListOrdered, referenceModel, comparisonName, analyses);
end

end


function project = runAllComparisons(project, modelListOrdered, referenceModel, comparisonName, analyses)
    % RUNALLCOMPARISONS Initializes the comparison slot and runs structural,
    % functional, and sampling comparisons as requested.
    % Structural comparison is always run (prerequisite for the others).
    
    % Initialisation
    project.comparisons.(comparisonName) = struct();
    project.comparisons.(comparisonName).modelNames = modelListOrdered;
    project.comparisons.(comparisonName).referenceModel = referenceModel;
    
    % Running structural comparison
    disp("Running structural comparison.");
    project.comparisons.(comparisonName).structuralComparison = structuralComparison(project, modelListOrdered, referenceModel);
    project.comparisons.(comparisonName).structuralAnalysisStatus = 1;
    
    % Functional comparison
    if any(matches(analyses, "functionalComparison"))
        disp("Running functional comparison.");
        project.comparisons.(comparisonName).functionalComparison.plots = functionalComparison(project, comparisonName);
    end
    
    % Sampling comparison
    if any(matches(analyses, "samplingComparison"))
        disp("Running sampling comparison.");
        project = samplingComparison(project, comparisonName);
    end

end