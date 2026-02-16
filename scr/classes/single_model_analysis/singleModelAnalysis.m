function project = singleModelAnalysis(project, modelList, analyses, parameterTable)
% this function runs the analysis on a model or a list or models and stores
% the results in a structure.
% Inputs: Project, Names of the models, List of wanted analysis (all by
% default)
% Available analysis:
%   FBA
%   FVA
%   shadow prices % need to be added
%   sampling
%   single_gene_deletion
%   double_gene_deletion
%   enrichment % need to be added
% Output : project with a field analysis

%% Check the structure of project
structureOk = checkStructureForSingleModelAnalysis(project, modelList);
if ~structureOk
    return
end
%% Check that given analysis are part of the list
%% Arguments (some are there by default)

%% Run analysis
for i = 1:numel(modelList)
    name = modelList{i};
    fprintf('Analysis running for model %s\n', name);
    
    if ~isfield(project.models.(name), 'analysis')
        project.models.(name).analysis = struct();
    end
    
    % Analysis id
    id = ['analysis_' char(datetime("now", "Format", "yyyyMMdd_HHmm"))];
    project.models.(name).analysis.(id) = struct();
    
    % Store the settings
    project.models.(name).analysis.(id).parameters = parameterTable;

    model = project.models.(name).model;
    
    % Searching for the objective function in the parameter table
    objFunction = parameterTable.Value(strcmp(parameterTable.Parameter, 'obj_function'));
    
    % Setting the objective function
    if isempty(objFunction)
        error('Parameter "obj_function" not found in the parameter table.');
        %return
    else
        model = changeObjective(model, objFunction);        
    end
    
    % FBA
    if any(strcmp(analyses, 'FBA'))
        disp('Running FBA.');
        
        % Loading parameters
        params = tableToParamsStruct(parameterTable, 'FBA', model);

        % Running FBA (specific format for arguments)
        FBA = optimizeCbModelParams(model, params);

        % Storing results
        project.models.(name).analysis.(id).FBA = FBA;
        %analysis.FBA = FBA;
    end

    % FVA
    if any(strcmp(analyses, 'FVA'))
        disp('Running FVA.');
        
        % Loading parameters
        params = tableToParamsStruct(parameterTable, 'FVA', model);
        
        % Running FVA
        [FVAmin, FVAmax] = fluxVariability(model, params);
        FVA = table(FVAmin, FVAmax, 'VariableNames', {'minFlux', 'maxFlux'});
        
        % Storing the results
        project.models.(name).analysis.(id).FVA = FVA;
        
    end
    
    % sampling
    if any(strcmp(analyses, 'sampling'))
        
        project.models.(name).analysis.(id).sampling = struct();
        
        % Loading parameters
        params = tableToParamsStruct(parameterTable, 'sampling', model);
        
        if ~isfield(params, 'osenseStr') && ~isfield(params, 'minNorm')
            if isfield(analysis, 'FBA')
                FBA = project.models.(name).analysis.(id).FBA;
            else
                FBA = optimizeCbModel(model, 'max', 'zero');
                project.models.(name).analysis.(id).FBA = FBA;
            end
        else
            FBA = optimizeCbModel(model, params.osenseStr, params.minNorm);
            project.models.(name).analysis.(id).FBA = FBA;
        end
        
        opt = FBA.f;
        disp("Setting the rxn bounds for the biomass_rxn to obj_threshold percent of the optimum")
        model = changeRxnBounds(model, objFunction, params.obj_threshold*opt, 'l'); % 0.9 usually
        
        if ~isfield(params, 'sampleFile') || (isfield(params, 'sampleFile') && isempty(params.sampleFile))
            sampleFile = char("sampleFile" + "_" + string(datetime("now", "Format", "yyyyMMdd_HHmm")));
        else
            sampleFile = char(string(params.path) + string(params.sampleFile) + "_" + string(datetime("now", "Format", "yyyyMMdd_HHmm")));
        end
        
        if ~isfield(params, 'samplerName') || (isfield(params, 'samplerName') && isempty(params.samplerName))
            samplerName = 'ACHR';
        else
            samplerName = params.samplerName;
            disp(samplerName);
        end
        
        if isfield(params, 'options')
            options = params.options;
            % Checking/setting required options
            if ~isfield(options, 'nPointsReturned')
                options.nPointsReturned = 2000;
            end
            if ~isfield(options, 'toRound') && samplerName == "CHRR"
                options.toRound = 1;
            end
            if ~isfield(options, 'nFiles')
                options.nFiles = 10;
            end
            if ~isfield(options, 'maxTime')
                options.maxTime = 36000;
            end
            if ~isfield(options, 'nWarmupPoints')
                options.nWarmupPoints = 2*size(model.S, 2);
            end
            if ~isfield(options, 'nStepsPerPoint')
                options.nStepsPerPoint = size(model.S, 2);
            end
        else
            options=[];
            options.nPointsReturned = 2000;
            options.nFiles = 10;  % increase this with the nPointsReturned (ratio 1 file ~ 100 samples)
            options.maxTime = 36000;  % 10 hours
            options.nWarmupPoints = 2*size(model.S, 2);
            options.nStepsPerPoint = size(model.S, 2);
            if samplerName == "CHRR"
                options.toRound = 1;
            end
        end
        
        
        disp("Running sampling.")
        if params.samplerName == "ACHR"
            [modelSampling, samples] = sampleCbModel(model, sampleFile, samplerName, options, model);
        elseif params.samplerName == "CHRR"
            [modelSampling, samples] = sampleCbModel(model, sampleFile, samplerName, options);
        else
            error("Currently this pipeline only allows the use of 2 samplers: ACHR and CHRR.")
        end
        disp("Sampling done.")
        
        save(params.path + "sampling_" + name + "_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".mat", 'modelSampling', 'samples');
        
        % Storing the results
        project.models.(name).analysis.(id).sampling.modelSampling = modelSampling;
        project.models.(name).analysis.(id).sampling.samples = samples;
        
    end
    
    % singleGeneDeletion
    if any(strcmp(analyses, 'singleGeneDeletion'))
        disp('Running singleGeneDeletion.');
        
        project.models.(name).analysis.(id).singleGeneDeletion = struct();
        
        % Loading parameters
        params = tableToParamsStruct(parameterTable, 'singleGeneDeletion', model);
        
        if ~isfield(params, 'method')
            method = 'FBA';
        else
            method = params.method;
        end
        
        if ~isfield(params, 'geneList')
            geneList = '';
            disp('All genes will be deleted by default.');
        else
            geneList = params.geneList;
        end
        
        if ~isfield(params, 'verbFlag')
            verbFlag = false;
        else
            verbFlag = params.verbFlag;
        end    
        
        if ~isfield(params, 'uniqueGene')
            uniqueGene = 0;
        else
            uniqueGene = params.uniqueGene;
        end
        
        [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model, method, geneList, verbFlag, uniqueGene);
        
        % Storing the results
        project.models.(name).analysis.(id).singleGeneDeletion.grRatio = grRatio;
        project.models.(name).analysis.(id).singleGeneDeletion.grRateKO = grRateKO;
        project.models.(name).analysis.(id).singleGeneDeletion.grRateWT = grRateWT;
        project.models.(name).analysis.(id).singleGeneDeletion.hasEffect = hasEffect;
        project.models.(name).analysis.(id).singleGeneDeletion.delRxns = delRxns;
        project.models.(name).analysis.(id).singleGeneDeletion.fluxSolution = fluxSolution;
        
    end
    
    % doubleGeneDeletion
    if any(strcmp(analyses, 'doubleGeneDeletion'))
        disp('Running doubleGeneDeletion.');
        
        project.models.(name).analysis.(id).doubleGeneDeletion = struct();
        
        % Loading parameters
        params = tableToParamsStruct(parameterTable, 'doubleGeneDeletion', model);
        
        if ~isfield(params, 'method')
            method = 'FBA';
        else
            method = params.method;
        end
        
        if ~isfield(params, 'geneList1') && ~isfield(params, 'geneList2')
            geneList1 = '';
            geneList2 = geneList1;
            disp('Every combination of genes will be deleted by default.');
        elseif isfield(params, 'geneList1') && ~isfield(params, 'geneList2')
            geneList1 = params.geneList1;
            geneList2 = geneList1;
            disp('By default, geneList2 will be equal to geneList1.');
        elseif ~isfield(params, 'geneList1') && isfield(params, 'geneList2')
            error("geneList1 is missing.")
        else
           geneList1 = params.geneList1;
           geneList2 = params.geneList2;
        end
       
        if ~isfield(params, 'printLevel')
            printLevel = 0;
        else
            printLevel = params.printLevel;
        end
        
        [grRatioDble, grRateKO, grRateWT] = doubleGeneDeletion(model, method, geneList1, geneList2, printLevel);
        
        % Storing the results
        project.models.(name).analysis.(id).doubleGeneDeletion.grRatioDble = grRatioDble;
        project.models.(name).analysis.(id).doubleGeneDeletion.grRateKO = grRateKO;
        project.models.(name).analysis.(id).doubleGeneDeletion.grRateWT = grRateWT;
        
    end
    
end

end
