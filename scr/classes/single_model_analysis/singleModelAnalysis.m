function project = singleModelAnalysis(project, modelList, analyses, parameterTable)
% this function runs the analysis on a model or a list or models and stores
% the results in a structure.
% Inputs: Project, Names of the models, List of wanted analysis (all by
% default)
% Available analysis:
%   FBA
%   FVA
%   sampling
%   single_gene_deletion
%   double_gene_deletion
%   enrichment ?
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
        project.models.(name).analysis.FBA = FBA;
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
        project.models.(name).analysis.FVA = FVA;
        
    end
    
    % sampling
    if any(strcmp(analyses, 'sampling'))
        
        project.models.(name).analysis.sampling = struct();
        
        % Loading parameters
        params = tableToParamsStruct(parameterTable, 'sampling', model);
        
        if ~isfield(params, 'osenseStr') && ~isfield(params, 'minNorm')
            if isfield(analysis, 'FBA')
                FBA = project.models.(name).analysis.FBA;
            else
                FBA = optimizeCbModel(model, 'max', 'zero');
            end
        else
            FBA = optimizeCbModel(model, params.osenseStr, params.minNorm);
        end
        
        opt = FBA.f;
        disp("Setting the rxn bounds for the biomass_rxn to obj_threshold percent of the optimum")
        model = changeRxnBounds(model, objFunction, params.obj_threshold*opt, 'l'); % 0.9 usually
        
        if ~isfield(params, 'sampleFile') || (isfield(params, 'sampleFile') && isempty(params.sampleFile))
            sampleFile = char("sampleFile" + "_" + string(datetime("now", "Format", "yyyyMMdd_HHmm")));
        else
            sampleFile = char(params.path + params.sampleFile + "_" + string(datetime("now", "Format", "yyyyMMdd_HHmm")));
        end
        
        if isfield(params, 'samplerName') || (isfield(params, 'samplerName') && isempty(params.samplerName))
            samplerName = 'ACHR';
        else
            samplerName = params.samplerName;
        end
        
        if isfield(params, 'options')
            options = params.options;
            % Checking/setting required options
            if ~isfield(options, 'nPointsReturned')
                options.nPointsReturned = 2000;
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
        end
        
        disp("Running sampling.")
        [modelSampling, samples] = sampleCbModel(model, sampleFile, samplerName, options, model);
        disp("Sampling done.")
        
        save(params.path + "sampling_" + name + "_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".mat", 'modelSampling', 'samples');
        
        % Storing the results
        project.models.(name).analysis.sampling.modelSampling = modelSampling;
        project.models.(name).analysis.sampling.samples = samples;
        
    end
    % singleGeneDeletion
    % doubleGeneDeletion
    
end

end
