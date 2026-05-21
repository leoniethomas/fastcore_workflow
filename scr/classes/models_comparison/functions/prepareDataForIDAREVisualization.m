function prepareDataForIDAREVisualization(project, comparisonName,folderPath,options)
    arguments
        project 
        comparisonName (1,1) string
        folderPath (1,1) string 
        options =[]
    end
    referenceModel = project.comparisons.(comparisonName).referenceModel;
    modelList = project.comparisons.(comparisonName).modelNames;
    
    % get unique identifier + create folder to store idare output in
    folderToStore = folderPath + filesep + datestr(now, 'yyyymmdd_HHMMSS');
    mkdir(folderToStore)
    storeModels = folderToStore + filesep + "models";
    storeData = folderToStore + filesep + "data";
    mkdir(storeModels)
    mkdir(storeData)
    

    % save all three models as xml files into the folder
    modelNames = project.comparisons.(comparisonName).modelNames;
    for modelIdx=1:length(modelNames)
        model = project.models.(modelNames(modelIdx)).model;
        modelFileName = storeModels + filesep + modelNames(modelIdx);
        modelFileNameOrig = storeModels + filesep + "orig_" +  modelNames(modelIdx);
        save(modelFileName + ".mat",'model');
        exportToXML(modelFileName + ".mat",modelFileName + ".xml",modelFileNameOrig + ".xml");
    end

    % save reference model
    model = project.models.(referenceModel).model;
    modelFileName = storeModels + filesep + referenceModel + "_reference";
    save(modelFileName + ".mat",'model');
    exportToXML(modelFileName + ".mat",modelFileName + ".xml",modelFileNameOrig + ".xml");

    % save the data that belongs to the models in the data folder, ready to
    % be load 

    % rxns data
    %------------------

    orderedMappingRxnMatrix = project.comparisons.(comparisonName).rxn_mapping_table;
    % fba + fbafluxsum 
    orderedFba = project.comparisons.(comparisonName).orderedFba;
    % fva + fvaSim
    replacementValue = "analysis.active.FVA.minFlux"; % get the fba solution values
    orderedFVAmin = getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacementValue);
    replacementValue = "analysis.active.FVA.maxFlux"; % get the fba solution values
    orderedFVAmax = getOrderedFeatureMatrix(project,modelList,"rxns",referenceModel,replacementValue);
    [fvaSim,orderedFvasim, ~] = computeFvaSimilarity(project,comparisonName);
    % sampling + sampling fluxsum 
    orderedSamples = project.comparisons.(comparisonName).orderedSamples;
    
    % build reaction dataset to load into cytoscape and use in IDARE
    prefix = 'rxn_idx_ind_model_';
    orderedMappingRxnMatrix.Properties.VariableNames = ...
        strcat(prefix, orderedMappingRxnMatrix.Properties.VariableNames);
    orderedPresenceRxns = orderedMappingRxnMatrix{:,:} >0;
    orderedOverallPresenceRxns = sum(orderedPresenceRxns,2);
    labels = project.comparisons.(comparisonName).sampleModelLabels;
    data   = orderedSamples;
    orderedMeanSampling =cell2mat(arrayfun(@(l) mean(data(:,labels==l),2), unique(labels, 'stable'), 'UniformOutput', false)); % unique uses alphabetic order for the output --> need to switch to unique(A, "stable") to preserve the order
    
    activeInFba = sum(orderedFba ~= 0,2);
    activeInSampling = sum(orderedSamples ~=0, 2);

    % combine into one dataframe 

    rxnData = [orderedFba,orderedPresenceRxns,orderedOverallPresenceRxns,activeInFba,activeInSampling,orderedFVAmin, orderedFVAmax,...
                orderedMeanSampling];
    rxnDataColNames = ["fba_" + modelList',"rxn_presence_" + modelList', "overall_rxns_presence","active_in_fba", "active_in_sampling", "FVA_min_" + modelList', "FVA_max_" + modelList',...
                          "mean_sampling_value_" + modelList'];
    rxnDataTable = array2table(rxnData, 'VariableNames',rxnDataColNames);

    rxnDataTable = [ rxnDataTable,orderedMappingRxnMatrix];
    rxnDataTable.rxnNameMatModel = rxnDataTable.Properties.RowNames;
    rxnDataTable.Properties.RowNames = strcat('R_', rxnDataTable.Properties.RowNames);
    rxnDataTable.Properties.RowNames = regexprep(rxnDataTable.Properties.RowNames, {'\[', '\]', '-'}, {'__91__','__93__','__45__'});
    
    rxnDataTable.isExchange = findExcRxns(project.models.(referenceModel).model);
    rxnDataTable.subsystem = string(project.models.(referenceModel).model.subSystems);
    rxnDataTable.symbolGprRules = string(cellfun(@(rxnName)getRxnSymbolRule(project.models.(referenceModel),...
                                                   rxnName),string(project.models.(referenceModel).model.rxns),'UniformOutput', false));
    
    
    rxnDataTable.RxnFormula = string(printRxnFormula(project.models.(referenceModel).model));

    rxnDataTable = addvars(rxnDataTable, string(rxnDataTable.Properties.RowNames), 'Before', 1, 'NewVariableNames', 'rxn_names');
    rxnDataTable.Properties.RowNames = {};  % remove old row names


    writetable(rxnDataTable,storeData + filesep + "reaction_data.xlsx");

    % metabolite data
    %------------------
    rxnCountPerMetConnectivity = sum(project.models.(referenceModel).model.S ~= 0, 2);
    orderedFbaFluxsum = getFluxsum(project,comparisonName,[],[],"orderedFba");
    orderedFbaFluxsum = orderedFbaFluxsum{:,:};
    orderedSamplesFluxsum = getFluxsum(project,comparisonName,[],[],"orderedSamples");
    orderedSamplesFluxsum = orderedSamplesFluxsum{:,:};
    orderedMeanSamplesFluxsum =cell2mat(arrayfun(@(l) mean(orderedSamplesFluxsum(:,labels==l),2), unique(labels), 'UniformOutput', false));
    
    orderedFbaFluxsum_presence = sum(orderedFbaFluxsum > 0,2);
    orderedSamplesFluxsumPresence = sum(orderedSamplesFluxsum >0,2);

    metNamesModel = project.models.(referenceModel).model.mets;
    metNamesCytoscape = strcat('M_', metNamesModel);
    metNamesCytoscape = regexprep(metNamesCytoscape, {'\[', '\]'}, {'__91__','__93__'});


    metTableColnames = ["met_names","met_names_matfile", "fba_fluxsum_presence", ...
                          "samples_fluxsum_presence", "orderedFba_fluxsum_" + modelList', ...
                          "mean_sample_fluxsum_" + modelList', "rxn_count_per_met_connectivity"];
    
    metDataTable = array2table([string(metNamesCytoscape) , string(metNamesModel), orderedFbaFluxsum_presence, ...
                                  orderedSamplesFluxsumPresence,orderedFbaFluxsum, ...
                                  orderedMeanSamplesFluxsum ,full(rxnCountPerMetConnectivity)],...
                                  "VariableNames",metTableColnames);

    writetable(metDataTable,storeData + filesep + "metabolite_data.xlsx");


    function exportToXML(matFile, xmlFile,xmlFile_orig)

        %pyenv('Version','/Users/leonie.thomas/miniconda3/envs/cobra_py/bin/python');
        pyenv('Version', '/Users/vanille.lejal/.conda/envs/python_v3.12/python.exe');
        py.importlib.import_module('cobra.io');
        
        model = py.cobra.io.load_matlab_model(matFile);
        py.cobra.io.write_sbml_model(model, xmlFile_orig);
    
    
        % Remove GPR rules from reactions
        rxns_iter = py.iter(model.reactions);
        while true
            try
                rxn = py.next(rxns_iter);
                rxn.gene_reaction_rule = '';  % remove GPR
            catch
                break;
            end
        end
    
        %% Remove compartment info from metabolites
        mets_iter = py.iter(model.metabolites);
        while true
            try
                met = py.next(mets_iter);
                met.compartment = '';             % remove compartment
            catch
                break;
            end
        end
    
        
        py.cobra.io.write_sbml_model(model, xmlFile);
        
        fprintf('Exported %s to %s\n', matFile, xmlFile);
    end


end