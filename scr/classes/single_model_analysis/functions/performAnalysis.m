function project = performAnalysis(project, parameterTable, modelName, toPerform, id)
%% Arguments block
arguments
    project struct
    parameterTable table
    modelName {mustBeText}
    toPerform cell {mustBeText}
    id {mustBeText}
end

%% Shortcut
model = project.models.(modelName).model;

%% Objective function 
% Searching for the objective function in the parameter table
objFunction = parameterTable.Value(strcmp(parameterTable.Parameter, 'objFunction')); % to remove

% Setting the objective function
if isempty(objFunction)
    error('Parameter "objFunction" not found in the parameter table.');
    %return
else
    model = changeObjective(model, objFunction);        
end

%% Perform Analysis
% FBA
if any(strcmp(toPerform, 'FBA'))
    disp('Running FBA.');
    
    % Loading parameters
    params = tableToParamsStruct(parameterTable, 'FBA', model);

    % Running FBA (specific format for arguments)
    FBA = optimizeCbModelParams(model, params);

    % Storing results
    project.models.(modelName).analysis.(id).FBA = FBA;
    %analysis.FBA = FBA;
end

% FVA
if any(strcmp(toPerform, 'FVA'))
    disp('Running FVA.');
    project.models.(modelName).analysis.(id).FVA = struct();
    
    % Loading parameters
    params = tableToParamsStruct(parameterTable, 'FVA', model);
    
    % Running FVA
    [FVAmin, FVAmax] = fluxVariability(model, params);
    FVA = table(FVAmin, FVAmax, 'VariableNames', {'minFlux', 'maxFlux'});
    
    % Storing the results
    project.models.(modelName).analysis.(id).FVA.minMaxFluxes = FVA;

    % check which rxns can loop without an input - tag the loops in the model
    modelTestLoops = changeRxnBounds(model, model.rxns(findExcRxns(model)), 0, 'b');
    [Vmin,Vmax] = fluxVariability(modelTestLoops);
    project.models.(modelName).analysis.(id).FVA.loopStatus = Vmin ~= 0 | Vmax ~= 0; % -> loops in the model
    
end

% sampling
if any(strcmp(toPerform, 'sampling'))
    
    project.models.(modelName).analysis.(id).sampling = struct();
    
    % Loading parameters
    params = tableToParamsStruct(parameterTable, 'sampling', model);
    mkdir(params.path);
    
    if ~isfield(params, 'osenseStr') && ~isfield(params, 'minNorm')
        if isfield(project.models.(modelName).analysis.(id), 'FBA')
            FBA = project.models.(modelName).analysis.(id).FBA;
        else
            FBA = optimizeCbModel(model, 'max', 'zero');
            project.models.(modelName).analysis.(id).FBA = FBA;
        end
    else
        FBA = optimizeCbModel(model, params.osenseStr, params.minNorm);
        project.models.(modelName).analysis.(id).FBA = FBA;
    end
    
    opt = FBA.f;
    disp("Setting the rxn bounds for the biomassRxn to objThreshold percent of the optimum")
    model = changeRxnBounds(model, objFunction, params.objThreshold*opt, 'l'); % 0.9 usually
    
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
        %changeCobraSolverParams('LP', 'feasTol', 1e-9);
        % we have the problem that the sampling does not investigate
        % the whole solution space, maybe too slow mixing ? too small
        % steps -> increase or decrease the threshold -> I think we
        % need to increase it, to allow faster mixing
        changeCobraSolverParams('LP','feasTol',1e-5);
        options.optPercentage = params.objThreshold*100;
        if any(strcmp(toPerform, 'FVA'))
            model.lb = FVA.minFlux; % constraining the sampling space by the FVA boundaries helps
            model.ub = FVA.maxFlux; % but that means that the threshold for the FVA is autoomatically applied to the sampling
        end
        options.nPointsReturned = round(options.nPointsReturned/options.countSampleProcesses);
        
        % start independend strains, in order to account for
        % randomness in the sampling 
        samples = [];  % vertical concatenation 
        for strain=1:options.countSampleProcesses
            disp("Running strain " + string(strain))
            rng(strain);  % set seed based on loop inde
            [modelSampling, samp] = sampleCbModel(model, sampleFile, samplerName, options);
            samples = [samples, samp];  % vertical concatenation
        end
    % 
    % elseif params.samplerName == "ADSB" | params.samplerName == "ll_ACHRB" | params.sampler == "EDHRB"
    %     option.numSamples = params.options.nPointsReturned;
    %     option.stepsPerPoint = params.options.nStepsPerPoint;
    %     %option.algorithm = params.samplerName;
    %     option.numDiscarded = params.options.nWarmupPoints;
    %     option.parallelFlag = true;
    %     samples = looplessFluxSampler(model, option);
    else
        error("Currently this pipeline does not allow the use of other samplers than ACHR and CHRR.")
    end
    disp("Sampling done.")
    
    save(params.path + "sampling_" + modelName + "_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".mat", 'modelSampling', 'samples');
    
    % Storing the results
    project.models.(modelName).analysis.(id).sampling.modelSampling = modelSampling;
    project.models.(modelName).analysis.(id).sampling.samples = single(samples);
end

% get rid of the thermodynamically infeasible loops using cycleFreeFlux
if any(strcmp(toPerform, 'loopless'))
    params = tableToParamsStruct(parameterTable, 'loopless', model);
    
    analysisIdSampling = params.samplingToUse;
            
    project.models.(modelName).analysis.(id).sampling = project.models.(modelName).analysis.(analysisIdSampling).sampling;
    samples = double(project.models.(modelName).analysis.(analysisIdSampling).sampling.samples);     

        
        samples = double(samples);
        evalc('initCobraToolbox()');% otherwise cycleFreeFlux function is not found
        step = 500;% RAM can handle up to 1000 per loop approx. with 24RAM machine
        n = size(samples, 2);
        Vthermo = zeros(size(samples));
        thermoFeas = zeros(size(samples));
        thermoStatusMatrix   = zeros(size(samples));  % default 0 = corrected
        neededAttempts   = zeros(size(samples,1),1);  % default 0 = corrected
        looplessStatus   = zeros(size(samples)); %  loopless reaction = 1, still looping = 0
        h = waitbar(0, 'Processing samples to get rid of thermodynamic infeasible loops...');counter = 0;N = numel(1:step:n);
        
        for idx = 1:step:n
            counter = counter + 1;
            cols  = idx:min(idx+step-1, n);
            chunk = samples(:, cols);
            % Vthermo gives you the corrected fluxes
            % thermoFeas gives you a boolean whether a flux was
            % already feasible 1 or it was not feasible -1
            % does not tell us which have been corrected, and whether
            % they are feasible now!!!
            evalc('[Vthermo(:, cols), thermoFeas(:, cols)] = cycleFreeFlux(chunk, repmat(model.c,1,numel(cols)), model)');
            % check how many loops are now in the solution 
            llChunk = Vthermo(:, cols);
            evalc('[~, looplessStatus(:, cols)] = cycleFreeFlux(llChunk, repmat(model.c,1,numel(cols)), model)');
            waitbar(counter / N, h);
        end  

        close(h);
        eta = getCobraSolverParams('LP', 'feasTol') * 10;
        fluxChangedBool    = abs(samples - Vthermo) >= eta;
        thermoStatusMatrix(logical(thermoFeas))                              =  1;   % already feasible
        thermoStatusMatrix(~logical(thermoFeas) & ~fluxChangedBool)          = -1;   % forced bounds

        project.models.(modelName).analysis.(id).sampling.cycleFreeFlux.samplesLl = single(Vthermo);
        project.models.(modelName).analysis.(id).sampling.cycleFreeFlux.thermoFeas = uint8(thermoFeas);
        project.models.(modelName).analysis.(id).sampling.cycleFreeFlux.sampleStatusAfterCorrection = uint8(thermoStatusMatrix);
        project.models.(modelName).analysis.(id).sampling.cycleFreeFlux.neededAttempts = uint8(neededAttempts);
        project.models.(modelName).analysis.(id).sampling.cycleFreeFlux.looplessStatus = uint8(looplessStatus);
    

    
end

% FDR correction for sampling results
if any(strcmp(toPerform, 'kld'))
    % Performing the FDR correction shown by Bruno G. Galuzzi et al 2024 
    %-> link to paper: https://www.sciencedirect.com/science/article/pii/S1532046424000157?via%3Dihub
    
    project.models.(modelName).analysis.(id).kld = struct();
    
    % Loading parameters
    params = tableToParamsStruct(parameterTable, 'kld', model);
    
    if any(strcmp(toPerform, 'sampling'))
        % add the sampling to be one of the sets 
        samplingMatrix = project.models.(modelName).analysis.(id).sampling.samples;
    else 
        samplingMatrix = [];
    end

    changeCobraSolverParams('LP','feasTol',1e-5);
    FBA = optimizeCbModel(model, "max");
    boundSampling = FBA.f * params.objThreshold;
    disp("Setting the rxn bounds for the biomassRxn to objThreshold percent of the optimum")
    model = changeRxnBounds(model, objFunction, boundSampling, 'l'); % 0.9 usually

    if any(strcmp(toPerform, 'FVA'))
        model.lb = FVA.minFlux; % constraining the sampling space by the FVA boundaries helps
        model.ub = FVA.maxFlux; % but that means that the threshold for the FVA is autoomatically applied to the sampling
    end
    
    [kldMatrix,...
        pValueKld,samplingSets,fdr] = performKdlDivergenceAnalysis(model,samplingMatrix,...
                                                 'nPointsReturned',params.nPointsReturned,'numberOfIndSamplings',params.numberOfIndSamplings);
    
    project.models.(modelName).analysis.(id).kld.samplingSets = samplingSets;
    project.models.(modelName).analysis.(id).kld.kldMatrix = kldMatrix;
    project.models.(modelName).analysis.(id).kld.pValueKld = pValueKld;
    project.models.(modelName).analysis.(id).kld.fdr = fdr;

end

% singleGeneDeletion
if any(strcmp(toPerform, 'singleGeneDeletion'))
    disp('Running singleGeneDeletion.');
    
    project.models.(modelName).analysis.(id).singleGeneDeletion = struct();
    
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
    project.models.(modelName).analysis.(id).singleGeneDeletion.grRatio = grRatio;
    project.models.(modelName).analysis.(id).singleGeneDeletion.grRateKO = grRateKO;
    project.models.(modelName).analysis.(id).singleGeneDeletion.grRateWT = grRateWT;
    project.models.(modelName).analysis.(id).singleGeneDeletion.hasEffect = hasEffect;
    project.models.(modelName).analysis.(id).singleGeneDeletion.delRxns = delRxns;
    project.models.(modelName).analysis.(id).singleGeneDeletion.fluxSolution = fluxSolution;
    
end

% doubleGeneDeletion
if any(strcmp(toPerform, 'doubleGeneDeletion'))
    disp('Running doubleGeneDeletion.');
    
    project.models.(modelName).analysis.(id).doubleGeneDeletion = struct();
    
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
    project.models.(modelName).analysis.(id).doubleGeneDeletion.grRatioDble = grRatioDble;
    project.models.(modelName).analysis.(id).doubleGeneDeletion.grRateKO = grRateKO;
    project.models.(modelName).analysis.(id).doubleGeneDeletion.grRateWT = grRateWT;
    
end

end