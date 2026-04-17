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

        % check which rxns can loop without an input - tag the loops in the model
        model_test_loops = changeRxnBounds(model, model.rxns(findExcRxns(model)), 0, 'b');
        [Vmin,Vmax] = fluxVariability(model_test_loops);
        project.models.(name).analysis.(id).loop_status = Vmin ~= 0 | Vmax ~= 0; % -> loops in the model
        
    end
    
    % sampling
    if any(strcmp(analyses, 'sampling'))
        
        project.models.(name).analysis.(id).sampling = struct();
        
        % Loading parameters
        params = tableToParamsStruct(parameterTable, 'sampling', model);
        mkdir(params.path);
        
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
            %changeCobraSolverParams('LP', 'feasTol', 1e-9);
            % we have the problem that the sampling does not investigate
            % the whole solution space, maybe too slow mixing ? too small
            % steps -> increase or decrease the threshold -> I think we
            % need to increase it, to allow faster mixing
            changeCobraSolverParams('LP','feasTol',1e-5);
            options.optPercentage = params.obj_threshold*100;
            if any(strcmp(analyses, 'FVA'))
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
            error("Currently this pipeline only allows the use of 2 samplers: ACHR and CHRR.")
        end
        disp("Sampling done.")
        
        save(params.path + "sampling_" + name + "_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".mat", 'modelSampling', 'samples');
        
        % Storing the results
        project.models.(name).analysis.(id).sampling.modelSampling = modelSampling;
        project.models.(name).analysis.(id).sampling.samples = samples;
        
        % get rid of the thermodynamically infeasible loops using cycleFreeFlux!
        if isfield(params, 'loopless') && params.loopless
            
            evalc('initCobraToolbox()');% otherwise cycleFreeFlux function is not found
            step = 1;% RAM can handle up to 1000 per loop approx. with 24RAM machine
            n = size(samples, 2);
            Vthermo_all = zeros(size(samples));
            thermo_feas_all = zeros(size(samples));
            thermoStatusMatrix   = zeros(size(samples));  % default 0 = corrected
            needed_attempts   = zeros(size(samples,1),1);  % default 0 = corrected
            loop_status   = ones(size(samples,1),1);  % default 1 = there are some loops left 
            %h = waitbar(0, 'Processing samples to get rid of thermodynamic infeasible loops...');counter = 0;N = numel(1:step:n);

            for idx = 1:step:n
                idx
                %counter = counter + 1;
                cols  = idx:min(idx+step-1, n);
                chunk = samples(:, cols);
                for_loop_break = 0;
                % Vthermo gives you the corrected fluxes
                % thermo_feas gives you a boolean whether a flux was
                % already feasible 1 or it was not feasible -1
                % does not tell us which have been corrected, and whether
                % they are feasible now!!!
                for attempt = 1:20 % sometimes the solver just crashes, so we retry a couple of times
                    [Vthermo_all(:,cols), thermo_feas_all(:,cols)] = cycleFreeFlux(chunk, repmat(model.c, 1, numel(cols)), model);
                    if any(Vthermo_all(:,cols))
                        needed_attempts(cols) = attempt;
                        for_loop_break = 1;
                        break; 
                    end% if any of the values are not 0 then the solver did not crash for those solutions -> so we continue with the next sample
                    
                end
                
                 % in the case where the for loop broke, a solution was found, so here we optimize this solution to have no infeasible loops 
                 % cause when running cycleFreeFlux - the optimization  sometimes leads to new infeasible loops being
                 % introduced, by looping over the solutions we try to get a solution which does not have any, or as least as possible loops
                 % in the publication they state the a solution is only  fluxfree if the solution we obtain will can be put into  cycleFreeflux 
                 % again and return the same output as input
                if for_loop_break == 1 && ~all(thermo_feas_all(:,cols)) && isfield(params, 'iterative_loop_correction') && params.iterative_loop_correction % only go into the loop if set so in the parameter file
                    for attempt = 1:100 % sometimes the solver just crashes, so we retry a couple of times
                        old_solution = Vthermo_all(:,cols);
                        old_feas = thermo_feas_all(:,cols);
                        [Vthermo_all(:,cols), thermo_feas_all(:,cols)] = cycleFreeFlux(Vthermo_all(:,cols), repmat(model.c, 1, numel(cols)), model);
                        attempt
                        sum(thermo_feas_all(:,cols),1)
                        if sum(thermo_feas_all(:,cols),1) == 0
                            Vthermo_all(:,cols) = old_solution;
                            thermo_feas_all(:,cols) = old_feas;
                        end
                        if all(thermo_feas_all(:,cols)) % if the found solution can be put into cyclefreeflux and get returned from the function unchanged then the solution is truly thermodynamically feasible 
                            loop_status(cols) = 0; % a loopless solution was found
                            break; 
                        end
                        
                    end
                end
                %waitbar(counter / N, h);
            end  
            
            close(h);
            eta = getCobraSolverParams('LP', 'feasTol') * 10;
            fluxChangedBool    = abs(samples - Vthermo_all) >= eta;
            thermoStatusMatrix(logical(thermo_feas_all))                              =  1;   % already feasible
            thermoStatusMatrix(~logical(thermo_feas_all) & ~fluxChangedBool)          = -1;   % forced bounds

            project.models.(name).analysis.(id).sampling.cycleFreeFlux.samples_ll = Vthermo_all;
            project.models.(name).analysis.(id).sampling.cycleFreeFlux.thermo_feas = thermo_feas_all;
            project.models.(name).analysis.(id).sampling.cycleFreeFlux.sample_status_after_correction = thermoStatusMatrix;
            project.models.(name).analysis.(id).sampling.cycleFreeFlux.needed_attempts = needed_attempts;
            project.models.(name).analysis.(id).sampling.cycleFreeFlux.loop_status = loop_status;
        end

        
    end

    % FDR correction for sampling results
    if any(strcmp(analyses, 'kld'))
        % Performing the FDR correction shown by Bruno G. Galuzzi et al 2024 
        %-> link to paper: https://www.sciencedirect.com/science/article/pii/S1532046424000157?via%3Dihub
        
        project.models.(name).analysis.(id).kld = struct();
        
        % Loading parameters
        params = tableToParamsStruct(parameterTable, 'kld', model);
        
        if any(strcmp(analyses, 'sampling'))
            % add the sampling to be one of the sets 
            sampling_matrix = project.models.(name).analysis.(id).sampling.samples;
        else 
            sampling_matrix = [];
        end
        [kld_matrix,...
            p_value_kld,sampling_sets,fdr] = perform_kdl_divergence_analysis(model,sampling_matrix,...
                                                     'nPointsReturned',params.nPointsReturned,'number_of_ind_samplings',params.number_of_ind_samplings);
        
        project.models.(name).analysis.(id).kld.sampling_sets = sampling_sets;
        project.models.(name).analysis.(id).kld.kld_matrix = kld_matrix;
        project.models.(name).analysis.(id).kld.p_value_kld = p_value_kld;
        project.models.(name).analysis.(id).kld.fdr = fdr;

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
