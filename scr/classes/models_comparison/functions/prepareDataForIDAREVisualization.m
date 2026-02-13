function prepareDataForIDAREVisualization(project, comparison_name,folder_path,options)
    arguments
        project 
        comparison_name (1,1) string
        folder_path (1,1) string 
        options =[]
    end
    reference_model = project.comparisons.(comparison_name).reference_model;
    modelList = project.comparisons.(comparison_name).modelNames;
    
    % get unique identifier + create folder to store idare output in
    folder_to_store = folder_path + filesep + datestr(now, 'yyyymmdd_HHMMSS');
    mkdir(folder_to_store)
    store_models = folder_to_store + filesep + "models";
    store_data = folder_to_store + filesep + "data";
    mkdir(store_models)
    mkdir(store_data)
    

    % save all three models as xml files into the folder
    model_names = project.comparisons.(comparison_name).modelNames;
    for model_idx=1:length(model_names)
        model = project.models.(model_names(model_idx)).model;
        model_file_name = store_models + filesep + model_names(model_idx);
        save(model_file_name + ".mat",'model');
        exportToXML(model_file_name + ".mat",model_file_name + ".xml");
    end

    % save reference model
    model = project.models.(reference_model).model;
    model_file_name = store_models + filesep + reference_model + "_reference";
    save(model_file_name + ".mat",'model');
    exportToXML(model_file_name + ".mat",model_file_name + ".xml");

    % save the data that belongs to the models in the data folder, ready to
    % be load 

    % rxns data
    %------------------

    ordered_mapping_rxn_matrix = project.comparisons.(comparison_name).rxn_mapping_table;
    % fba + fbafluxsum 
    ordered_fba = project.comparisons.(comparison_name).ordered_fba;
    % fva + fva_sim
    replacement_value = "analysis.FVA.minFlux"; % get the fba solution values
    ordered_FVAmin = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
    replacement_value = "analysis.FVA.maxFlux"; % get the fba solution values
    ordered_FVAmax = getOrderedFeatureMatrix(project,modelList,"rxns",reference_model,replacement_value);
    [fva_sim,ordered_fvasim, ~] = compute_fva_similariy(project,comparison_name);
    % sampling + sampling fluxsum 
    ordered_samples = project.comparisons.(comparison_name).ordered_samples;
    
    % build reaction dataset to load into cytoscape and use in IDARE
    prefix = 'rxn_idx_ind_model_';
    ordered_mapping_rxn_matrix.Properties.VariableNames = ...
        strcat(prefix, ordered_mapping_rxn_matrix.Properties.VariableNames);
    ordered_presence_rxns = ordered_mapping_rxn_matrix{:,:} >0;
    ordered_overall_presence_rxns = sum(ordered_presence_rxns,2);
    labels = project.comparisons.(comparison_name).sample_model_labels;
    data   = ordered_samples;
    ordered_mean_sampling =cell2mat(arrayfun(@(l) mean(data(:,labels==l),2), unique(labels), 'UniformOutput', false));
    
    active_in_fba = sum(ordered_fba ~= 0,2);
    active_in_sampling = sum(ordered_samples ~=0, 2);

    % combine into one dataframe 

    rxn_data = [ordered_fba,ordered_presence_rxns,ordered_overall_presence_rxns,active_in_fba,active_in_sampling,ordered_FVAmin, ordered_FVAmax,...
                ordered_mean_sampling];
    rxn_data_col_names = ["fba_" + modelList',"rxn_presence_" + modelList', "overall_rxns_presence","active_in_fba", "active_in_sampling", "FVA_min_" + modelList', "FVA_max_" + modelList',...
                          "mean_sampling_value_" + modelList'];
    rxn_data_table = array2table(rxn_data, 'VariableNames',rxn_data_col_names);

    rxn_data_table = [ rxn_data_table,ordered_mapping_rxn_matrix];
    rxn_data_table.rxn_name_mat_model = rxn_data_table.Properties.RowNames;
    rxn_data_table.Properties.RowNames = strcat('R_', rxn_data_table.Properties.RowNames);
    rxn_data_table.Properties.RowNames = regexprep(rxn_data_table.Properties.RowNames, {'\[', '\]', '-'}, {'__91__','__93__','__45__'});
    
    rxn_data_table.is_exchange = findExcRxns(project.models.(reference_model).model);
    rxn_data_table.subsystem = string(project.models.(reference_model).model.subSystems);
    rxn_data_table.symbol_gpr_rules = string(cellfun(@(rxnName)get_rxn_symbol_rule(project.models.(reference_model),...
                                                   rxnName),string(project.models.(reference_model).model.rxns),'UniformOutput', false));
    
    
    rxn_data_table.RxnFormula = string(printRxnFormula(project.models.(reference_model).model));

    rxn_data_table = addvars(rxn_data_table, string(rxn_data_table.Properties.RowNames), 'Before', 1, 'NewVariableNames', 'rxn_names');
    rxn_data_table.Properties.RowNames = {};  % remove old row names


    writetable(rxn_data_table,store_data + filesep + "reaction_data.xlsx");

    % metabolite data
    %------------------
    rxn_count_per_met_connectivity = sum(project.models.(reference_model).model.S ~= 0, 2);
    ordered_fba_fluxsum = get_fluxsum(project,comparison_name,[],[],"ordered_fba");
    ordered_fba_fluxsum = ordered_fba_fluxsum{:,:};
    ordered_samples_fluxsum = get_fluxsum(project,comparison_name,[],[],"ordered_samples");
    ordered_samples_fluxsum = ordered_samples_fluxsum{:,:};
    ordered_mean_samples_fluxsum =cell2mat(arrayfun(@(l) mean(ordered_samples_fluxsum(:,labels==l),2), unique(labels), 'UniformOutput', false));
    
    ordered_fba_fluxsum_presence = sum(ordered_fba_fluxsum > 0,2);
    ordered_samples_fluxsum_presence = sum(ordered_samples_fluxsum >0,2);

    met_names_model = project.models.(reference_model).model.mets;
    met_names_cytoscape = strcat('M_', met_names_model);
    met_names_cytoscape = regexprep(met_names_cytoscape, {'\[', '\]'}, {'__91__','__93__'});


    met_table_colnames = ["met_names","met_names_matfile", "fba_fluxsum_presence", ...
                          "samples_fluxsum_presence", "ordered_fba_fluxsum_" + modelList', ...
                          "mean_sample_fluxsum_" + modelList', "rxn_count_per_met_connectivity"];
    
    met_data_table = array2table([string(met_names_cytoscape) , string(met_names_model), ordered_fba_fluxsum_presence, ...
                                  ordered_samples_fluxsum_presence,ordered_fba_fluxsum, ...
                                  ordered_mean_samples_fluxsum ,full(rxn_count_per_met_connectivity)],...
                                  "VariableNames",met_table_colnames);

    writetable(met_data_table,store_data + filesep + "metabolite_data.xlsx");

end

function exportToXML(matFile, xmlFile)

    pyenv('Version','/Users/leonie.thomas/miniconda3/envs/cobra_py/bin/python');
    py.importlib.import_module('cobra.io');
    
    model = py.cobra.io.load_matlab_model(matFile);


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








