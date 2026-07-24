function [project, comparisonName] = modelsComparison(project, modelList, analysisID, referenceModel, analyses, identifier)
    % This function runs a set of analysis for the comparison of the
    % specified models.
    % A number of analysis are run: 
    % - structural analysis: based on the differential presence of
    %   metabolites, genes and reactions in the different models 
    % - functional analysis: based on the quantitative values like FVA,
    %   FBA
    % - sampling analysis: comparative analysis based on the sampling
    %   results
    % 
    % Inputs: 
    %   - project:          the object which is the output of the single_model_analysis
    %                       entailing the results of fba,fva,sampling, single gene
    %                       deletion etc. for a single model 
    %   - modelList:        the list of Model names to be included in the comparison 
    %   - analysisID:       the individual analysis slots used for the
    %                       comparison (as a default the analysis performed
    %                       most recent is choosen)
    %   - referenceModel:  the reference model used to compute the relative reaction presence
    %   - analyses:         the list of analyses which should be performed 
    %                       + modelStructureComparison: investigates the differences between the models 
    %                         on a structural level, is a gene,metabolite,or rxns present or not ? 
    %                       + modelFunctionalComparison: investigates the differences between models on a functional
    %                         level, how much flux do reactions carry in FVA, FBA solutions? 
    %                       + samplingComparison: investigates the differences between models sampling solutions,
    %                         investigates samples solution space 
    %                         
    %   - identifier:       a string, will be added as a postfix to the analysis name, can be choosen freely
    %                       default)
    %   
    % Output : 
    %   - project:          project object with a added comparison field entailing
    %                       all the output, modelcomparison information
    %   - comparisonName:   gives back the name of the comparison name added
    arguments
        project (1,1) struct 
        modelList (1,:) string
        analysisID (1,:) string
        referenceModel (1,1) string
        analyses  (1,:) string {mustBeMember(analyses, {'structuralComparison', 'functionalComparison', 'samplingComparison' , 'IDAREoutput'})} = 'structural'
        identifier (1,1) string = string(datetime('now', 'Format', '_yyyyMMdd_HHmmss'))
    end
    
    %% Check project and models format for comparison

    % Check first the reference model
    checkProjectFormat(project, referenceModel)

    % Check then models to compare
    checkProjectFormat(project, modelList, repmat({'active'}, 1, numel(modelList)))

    %% Create comparison slot
    % give the comparison the name of all models + a identifier choosen

    order = string(fieldnames(project.models));
    [~, idx] = ismember(modelList, order);
    [~, sortIdx] = sort(idx);
    
    modelListOrdered = modelList(sortIdx);
    clear modelList
    comparisonName = join(modelListOrdered, "_vs_") + "__" + identifier;
    % does this comparison object already exist, was the structural
    % analysis performed ? The structural analysis is needed to do the two
    % other analysis, therefore we check if it was already run, to not run
    % it again if it was already created
    if isstruct(project.comparisons)
        fields = fieldnames(project.comparisons);
    else
        disp('project.comparisons is empty or not a struct');
        fields = {};
    end

    if ismember(comparisonName,fields)
        structureAnalysisAlreadyRun = project.comparisons.(comparisonName).structure_analysis_status;
    else
        structureAnalysisAlreadyRun = 0;
    end
    
    % run structure comparison if: it was the input of the function, if
    % another comparison was input but the structure analysis was not run
    % yet (no comparison object in the comparisons slot, or a comparison
    % object that is not complete -> structure_analysis_already_ran  not
    % defined
    if any(matches(analyses, "modelStructureComparison")) || ~structureAnalysisAlreadyRun || ~ismember(comparisonName,string(fieldnames(project.comparisons)))
        project.comparisons.(comparisonName).modelList = modelListOrdered;
        project.comparisons.(comparisonName).referenceModel = referenceModel;
        project.comparisons.(comparisonName).analysisID = analysisID;
        
        % -- run structural comparison - always has to be run 
    
        project.comparisons.(comparisonName) = modelStructuralComparison(project,modelListOrdered,referenceModel);
        project.comparisons.(comparisonName).referenceModel = referenceModel;
        project.comparisons.(comparisonName).structure_analysis_status = 1;
    end

    % -- run functional comparions
    if any(matches(analyses, "modelFunctionalComparison"))
        project.comparisons.(comparisonName).plots.funct = modelFunctionalComparison(project, comparisonName);
    end

    % -- run sampling comparison
    if any(matches(analyses, "modelSamplingComparison"))
        project = modelSamplingComparison(project,comparisonName);
    end

end





