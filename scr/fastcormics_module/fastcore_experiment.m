classdef fastcore_experiment
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        orig_model
        consistent_medium_constrained_model
        condition_specific_models
        dico
        medium
        data
        optional_settings
        script_parameters
    end
    
    methods
        %%
        
        function obj = fastcore_experiment(def_run_file)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.script_parameters = read_in_run_def_file(def_run_file);


        end

        function obj = add_models(obj,model,dico)
            disp("Running fastcc and getting rid of unconsistent rxns!")
            A = fastcc(model, 1e-4, 1);
            % 
            % % remove non consistent reactions from model
            model=removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));
            % % check if the biomass reactions are still there
            consistency_check(model);
            
            obj.orig_model = model;
            obj.dico = dico;
            
        end
        function obj = add_expression_data(obj,data, disc_data)
            obj.data = data;
            obj.data.source = disc_data;
        end
        
        
        function obj = medium_constrain(obj,mode,close_all_exchange_rxns, column_media_rxn_abbr)
            
            arguments
                obj
                mode (1,1) string {mustBeMember(mode,["set_fluxes","set_concentrations"])} ="set_fluxes"
                close_all_exchange_rxns (1,1) double {mustBeMember(close_all_exchange_rxns,[0 1])} =0
                column_media_rxn_abbr (1,1) string ="ExRxns_Recon3D"
            end
            medium_data = obj.medium;
            
            % find rxns defined in media composition
            model = obj.orig_model;
            
            [~,idx, idx_fluxes_in_model] = intersect(medium_data.medium_composition.(column_media_rxn_abbr),...
                                                     model.rxns); 
                                                 
            % set fluxes for metabolites in the medium
            if mode == "set_fluxes"
                model.lb(idx_fluxes_in_model) = medium_data.medium_composition{idx,"Flux_mmol_gDW_h"};
            elseif mode == "set_concentrations"
                model.lb(idx_fluxes_in_model) = medium_data.medium_composition{idx,"Concentration_M"};
            end
            
            [EX, UPT] = findExcRxns(model);
            needed_mets = medium_data.manual_set_boundaries.wanted_import;
            Ex_to_close = setdiff(model.rxns(findExcRxns(model)),...
                                             [medium_data.medium_composition.(column_media_rxn_abbr);...
                                              findRxnsFromMets(model, needed_mets)]);

            % set the rxns lower and upper bound, rxns that we set that we do not want
            % to have, are those also reasonable in my case, for my data ? 
            model.ub(find(ismember(model.rxns,split(medium_data.manual_set_boundaries.unwanted_export, ";"))))=0; 
            model.lb(find(ismember(model.rxns,split(medium_data.manual_set_boundaries.unwanted_import, ";"))))=0; 

            
            if close_all_exchange_rxns
                disp("All the exchange rxns not defined in the medium will be closed!")
                model.lb(findRxnIDs(model, Ex_to_close))=0; 
            end 
            
            
            
            disp("Running fastcc and getting rid of unconsistent rxns!")
            A = fastcc(model, 1e-4, 1);

            % remove non consistent reactions from model
            model=removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));
            % check if the biomass reactions are still there
            consistency_check(model);
            
            obj.consistent_medium_constrained_model = model;
        end

         function obj = add_sampling_to_fastcore_experiment(obj,sampling_files,run_fluxsum)
             %UNTITLED4 Construct an instance of this class
             %   Detailed explanation goes here
             arguments
             obj
             sampling_files
             run_fluxsum =0
             end
             sample_data = cell2struct(arrayfun(@(x) load(fullfile(x)), sampling_files','UniformOutput', false),...
                           regexprep(sampling_files, ".mat", ""),...
                           2);
             models =  structfun(@(y) y.x.modelSampling,sample_data,'UniformOutput', false);
             samples =  structfun(@(y) y.x.samples,sample_data,'UniformOutput', false);
             
             sample_names = arrayfun(@(x) repmat(x, 1,size(samples.(string(x)),2)) , fieldnames(samples)'  , 'UniformOutput' , false);
             obj.sample_labels = string([sample_names{:}]);
             
             obj.fastcore_runs = cell2struct(arrayfun(@(x) fastcore_run(models.(x),samples.(x),run_fluxsum),...
                                             string(fieldnames(models)),'UniformOutput', false),...
                                             regexprep(sampling_files, ".mat", ""),...
                                             1);
                                         
             obj.run_names = string(fieldnames(obj.fastcore_runs));
             
         end
         
         function obj = add_fba_to_experiment(obj)
             
         end
        
        function obj = join_sampling_output(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            all_rxns = cellfun(@(x) obj.fastcore_runs.(x).model.rxns, string(fieldnames(obj.fastcore_runs)),'UniformOutput',false);
            all_rxns = unique(vertcat(all_rxns{:}));

            samples_ordered = arrayfun(@(x) get_sampling_orig_order(obj.fastcore_runs.(x).model,obj.fastcore_runs.(x).sampling,all_rxns), ...
                                              string(fieldnames(obj.fastcore_runs)),...
                                              'UniformOutput',false);
            biomass_idx = find(ismember(all_rxns, "biomass_reaction"))

            obj.samples = cell2mat(samples_ordered');
            obj.rxn_names = all_rxns;
            
            disp("Biomass rxn value min and max overall samples:")
            min(obj.samples(biomass_idx,:))
            max(obj.samples(biomass_idx,:))

        end
        function [min_out,max_out] = join_fva_output(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            %all_rxns = cellfun(@(x) obj.condition_models.(x).rxns, string(fieldnames(obj.condition_models)),'UniformOutput',false);
            %all_rxns = unique(vertcat(all_rxns{:}));
            all_rxns = obj.orig_model.rxns;

            max_ordered = arrayfun(@(x) get_sampling_orig_order(obj.condition_models.(x),obj.fva.maxFlux.(x),all_rxns), ...
                                              string(fieldnames(obj.condition_models)),...
                                              'UniformOutput',false);
            max_out = array2table(cell2mat(max_ordered'), 'VariableNames', fieldnames(obj.condition_models));
            max_out.Properties.RowNames = all_rxns;
            min_ordered = arrayfun(@(x) get_sampling_orig_order(obj.condition_models.(x),obj.fva.minFlux.(x),all_rxns), ...
                                              string(fieldnames(obj.condition_models)),...
                                              'UniformOutput',false);

            min_out = array2table(cell2mat(min_ordered'), 'VariableNames', fieldnames(obj.condition_models));
            min_out.Properties.RowNames = all_rxns;
        end

        function obj = join_fluxsum_output(obj)
            %all_mets = cellfun(@(x) obj.fastcore_runs.(x).model.mets, string(fieldnames(obj.fastcore_runs)),'UniformOutput',false);
            %all_mets = unique(vertcat(all_mets{:}));
            all_mets = obj.orig_model.mets;
            
            %x = 'samplingResults_MDA_MB231_Cont_NO_model_20250602_090252';
            fluxsum = arrayfun(@(x) get_sampling_orig_order_mets(obj.fastcore_runs.(x).model,...
                                                                         obj.fastcore_runs.(x).fluxsum,...
                                                                         all_mets),...
                                       string(fieldnames(obj.fastcore_runs)),...
                                       'UniformOutput',false);
                                   
            obj.fluxsum = cell2mat(fluxsum');
            obj.met_names = all_mets;
            
        end
        function obj = change_model_labels(obj,keys,values)

            for i = 1:numel(keys)
                obj.run_names(obj.run_names == keys(i)) = values(i);
                obj.sample_labels(obj.sample_labels == keys(i)) = values(i);
            end
        end
        %%        
        function [obj,mean_sil,homogen] = visualize_sampling(obj, num_clusters_kmeans, ...
                                                                pc_x, pc_y,disp_fig,data_slot_to_use,use_pca_slot)
                                                            
            arguments
                obj
                num_clusters_kmeans
                pc_x
                pc_y
                disp_fig
                data_slot_to_use ="samples"
                use_pca_slot =0
            end
                                                            

            samples = obj.(data_slot_to_use)';
            
            % check if the pca was already performed!
            if ~use_pca_slot
                disp("compute pcs!")
                [coeff,score,latent,tsquared,explained] = pca(samples);
                obj.transformed_samples.pca.score = score;
                obj.transformed_samples.pca.explained = explained;
           else
                disp("using the pcs which were already computed!")
                score = obj.transformed_samples.pca.score;
                explained = obj.transformed_samples.pca.explained;
            end
            
            if num_clusters_kmeans
                disp("compute kmeans!")
                km = kmeans(samples,num_clusters_kmeans);
                [s,~] = silhouette(samples,km,'Euclidean');
                mean_sil = mean(s);
                %%
                homogen = [];
                for k = unique(km)'
                    labels_in_cluster = obj.sample_labels(find(km == k));
                    [s,~,j]=unique(labels_in_cluster);
                    f = s{mode(j)};
                    m = sum(labels_in_cluster == f);
                    homogen = [ homogen, m/length(labels_in_cluster)];
                end
                homogen = mean(homogen);
                
                if disp_fig
                    figure
                    scatter(score(:,pc_x),score(:,pc_y),5,km)

                    for i = unique(km)'
                        disp(i)
                        idx_label_samples = find(km == i);
                        %disp(mean(idx_label_samples))
                        pos_text = mean(score(idx_label_samples,:));
                        text(pos_text(pc_x),pos_text(pc_y), num2str(i), 'VerticalAlignment',...
                             'bottom', 'HorizontalAlignment', 'right');
                        hold on 
                    end

                    title("PC "+ num2str(pc_x) +  " & " +  num2str(pc_y) + " with the cluster ids from kmeans - " + data_slot_to_use)
                    xlabel("PC" + num2str(pc_x) + " var: " + explained(pc_x))
                    ylabel("PC" + num2str(pc_y) + " var: " + explained(pc_y))
                    hold off
                end
                
                figure
                hist(categorical(string(obj.sample_labels)' + "__" + km))
                title("Count of samples per condition - kmeans cluster id combination - " + data_slot_to_use)
                xlabel("condition__cluster_id")
                ylabel("count of samples")
                xtickangle(45)  
            end
            
            if disp_fig
                figure
                %id_biomass_ordered = 3625;
                %scatter(score(:,pc_x),score(:,pc_y),5,samples(:,id_biomass_ordered))
                scatter(score(:,pc_x),score(:,pc_y),5,categorical(obj.sample_labels))

                unique_labels = categories(categorical(obj.sample_labels));
                for i = 1:numel(unique_labels)

                    l = string(unique_labels(i));
                    idx_label_samples = find(obj.sample_labels == l);
                    %disp(mean(idx_label_samples))
                    pos_text = mean(score(idx_label_samples,:));
                    text(pos_text(pc_x),pos_text(pc_y), unique_labels(i), ...
                         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
                    hold on 
                end
                title("PC "+ num2str(pc_x) +  " & " +  num2str(pc_y) + " with the condition labels - " + data_slot_to_use)
                xlabel("PC" + num2str(pc_x) + " var: " + explained(pc_x))
                ylabel("PC" + num2str(pc_y) + " var: " + explained(pc_y))
                hold off
                
            end

            %%
            
        end 
        
        
        function [] = visualize_sampling_rxn_distribution(obj,rxn_id,data_slot_to_use)
            
            arguments
                obj
                rxn_id
                data_slot_to_use ="samples" 
            end
            
            data = obj.(data_slot_to_use);
            rxn_samp_fluxes = obj.samples(rxn_id,:);
            plotting_range_min = min(data(rxn_id,:));
            plotting_range_max = max(data(rxn_id,:));
            plotting_range = linspace(plotting_range_min,plotting_range_max,50);

            figure
            for i=unique(obj.sample_labels)
                idx = find(obj.sample_labels == i);
                %[probability_estimate,xi] = ksdensity(rxn_samp_fluxes(1,idx));
                %plot(xi,probability_estimate*100,'LineWidth',1); % multiplied with 100 to have %
                %trapz(probability_estimate, xi) % should approx to 1 % integral should sum up to 1
                [y,x] = hist(rxn_samp_fluxes(1,idx),plotting_range)
                plot(x,(y/length(idx))*100,'LineWidth',1);
                hold on
            end
            legend(unique(obj.sample_labels))
            xlabel("rxn flux value")
            ylabel("probability of obtaining x [%]")
            title("Probability distributions between different models given the performed sampling - rxn: " +  obj.rxn_names(rxn_id) + " (idx: " + num2str(rxn_id) + " )" )
            hold off
            
        end
        
        function obj = compute_fva_similariy(obj)
            modelNames = fieldnames(obj.condition_models);
            n = numel(modelNames);
            
            obj.fva_similarity = eye(n); % Diagonal is 1
            obj.fva_similarity_rxns = cell(n);
            
            fvaData = cell(1, n);
            for i = 1:n
                fvaData{i} = [obj.fva.minFlux{:, i}, obj.fva.maxFlux{:, i}];
            end

            
            indexPairs = find(triu(true(n), 1)); % Upper triangle linear indices
            [rowIdx, colIdx] = ind2sub([n, n], indexPairs);

            % Step 6: Functional-style loop without nested for-loops
            for k = 1:length(indexPairs)
                y = rowIdx(k);
                x = colIdx(k);

                [overallSim, rxnSim] = FVAsimilarity(fvaData{y}, fvaData{x});

                % Fill both [y,x] and [x,y] due to symmetry
                obj.fva_similarity(y,x) = overallSim;
                obj.fva_similarity(x,y) = overallSim;

                obj.fva_similarity_rxns{y,x} = rxnSim;
                obj.fva_similarity_rxns{x,y} = rxnSim;
            end
        end

        
        function stats = diff_flux_testing(obj,models_to_compare,writetofile,figflag,filepath,filename)
            arguments
                obj
                models_to_compare (1,2) string {mustBeMember_inModel(models_to_compare,obj)}
                writetofile (1,1) double {mustBeMember(writetofile,[1,0])} =1
                figflag (1,1) double {mustBeMember(figflag,[1,0])} =1
                filepath (1,1) string ="./"
                filename (1,1) string {mustBeFileType(filename,"xlsx")} =strjoin(["stats", models_to_compare,".xlsx"],"_")  
            end
            
            
            model1 = models_to_compare(1);
            model2 = models_to_compare(2);
            samples_model1 = obj.samples(:,obj.sample_labels == model1);
            samples_model2 = obj.samples(:,obj.sample_labels == model2);


            stats = zeros(size(obj.samples, 1), 1);  % Preallocate for speed
            for counter=1:size(obj.samples, 1)
                rxn_sample_value_model1 = samples_model1(counter,:);
                rxn_sample_value_model2 = samples_model2(counter,:);    
                stats(counter) = ranksum(rxn_sample_value_model1,rxn_sample_value_model2 );
            end

            % adjust p-value for multiple testing
            stats = [stats,mafdr(stats)];

            % compute signal to noise ration 
            % + adding 1000 to all samples flux values, so that we are only dealing on
            % a positive scale 
            % + we are doing the snr ratio -> deviding through the std of both
            % distributions -> normally only the noise std 
            % + but here we do not really have a baseline and a noise distribution 
            mean_model1 = mean(samples_model1 + 1000, 2)-1000;  % mean of each row (observation)
            median_model1 = median(samples_model1 + 1000, 2)-1000;  % mean of each row (observation)
            std_model1 = std(samples_model1 + 1000, 0, 2);  % std of each row (observation)

            mean_model2 = mean(samples_model2 + 1000, 2)-1000;  % mean of each row (observation)
            median_model2 = median(samples_model2 + 1000, 2)-1000;  % mean of each row (observation)
            std_model2 = std(samples_model2 + 1000, 0, 2);  % std of each row (observation)

            snr = (mean_model2 - mean_model1) ./ (std_model1 + std_model2);

            log2FC=log2(abs(mean_model2./mean_model1));
            stats =[mean_model1, mean_model2,median_model1, median_model2, std_model1, std_model2, ...
                          mean_model1 - mean_model2, log2FC, snr, ...
                          stats];
            
            if writetofile
                stats_tab=array2table(stats,'RowNames',obj.rxn_names,'VariableNames',{'mean_model1','mean_model2','median_model1','median_model2','std_model1','std_model2','diff','log2FC','SNR','pValue','p_adj'});
                writetable(stats_tab, filepath + filename,'WriteRowNames',true)
            end
            if figflag
                figure
                histogram(-log2(stats(:,11)))
                %title('P values (-log10)')
                title( 'Histogram P adj values (-log2)')
                
                figure
                %hist(log10(stats(:,4)))
                histogram(stats(:,8))
                %title('log10 foldchange (mean(B)/mean(A))')
                title('Histogram log2 foldchange mean\_model2/mean\_model1))')

                figure
                % plot(log10(stats(:,4)),stats(:,6),'*')
                % title('vulcano: log10 foldchange vs -log10(P)')
                plot(stats(:,8),-log2(stats(:,11)),'.')
                hold on 
                yline(-log2(0.001))
                xline(1.5)
                xline(-1.5)
                text(-10,100,["adj\_pval < 0.001", "abs(logFC) > 1.5"])
                title('vulcano: log2 foldchange vs -log2(P\_adj)')
                xlabel("log2FC flux")
                ylabel("-log2(p\_adj)")
            end

        end
        
        
    %%%%%%%%%%%%% medium functions!!

        function obj = add_medium(obj,filenames,filelocation)
            %arguments
            %       filenames (1,:) string {mustBeFileType(filenames,"tsv")}
            %       filelocation (1,:) string
            %end
            obj.medium.file_names = filenames;
            
            if length(filenames) == length(filelocation)
                obj.medium.file_location = filelocation;
            elseif length(filelocation) ==1
                obj.medium.file_location = filelocation(ones(1,length(filenames)));
            else
                error("The filelocation arguments needs to be one string, or as many strings as files are given!")
            end     
        end
        
        function obj = read_medium_files(obj, column_to_join_by,columns_in_both)
           arguments
                 obj
                 column_to_join_by
                 columns_in_both =["Mets_Recon3D","ExRxns_Recon3D","ExRxns_HumanGEM","Mets_HumanGEM"]
           end
           col_postfix = regexprep(obj.medium.file_names,".tsv","");
           column_to_join = column_to_join_by + "_" + col_postfix;
           med = {};
           for idx_file=1:length(obj.medium.file_names)
               file = obj.medium.file_location(idx_file) + filesep + obj.medium.file_names(idx_file);
               med{idx_file} = readtable(file, ...
                                         'FileType', 'text', 'Delimiter', '\t');               
               med{idx_file}.source = file(ones(size(med{idx_file},1),1))
               
               med{idx_file}.Properties.VariableNames(string(med{idx_file}.Properties.VariableNames) == column_to_join_by(idx_file)) = column_to_join(idx_file);
           end
           
           tab = med{1};
           col_conc = column_to_join(1);
           for med_idx=2:length(obj.medium.file_names)
                medium = outerjoin(tab,med{med_idx},...
                                   'MergeKeys',true,'Keys',...
                                   columns_in_both);
                                    
                medium{find(isnan(medium{:, [col_conc]})),col_conc} = 0;
                medium{find(isnan(medium{:, [column_to_join(med_idx)]})),column_to_join(med_idx)} = 0;
                
                % adding up the metabolites, in case they are in both media components
                medium{:, ["Concentration_M"]} = medium{:, [col_conc]} + medium{:, [column_to_join(med_idx)]};
                medium{:, ["Concentration_mM"]} = medium{:, ["Concentration_M"]}*1000;
                medium.(col_conc) = [];
                medium.(column_to_join(med_idx)) = []; 
                col_conc = "Concentration_M_tab";
                medium.(col_conc) = medium.Concentration_M;
                tab = medium;
           end
           medium.(col_conc) = [];
           obj.medium.medium_composition = medium;
            
        end
        
        function obj = add_additional_rxns_boundaries(obj, unwanted_import, unwanted_export,wanted_import, wanted_export)
           arguments
              obj
              unwanted_import (1,:) string =[]
              unwanted_export (1,:) string =[]
              wanted_import (1,:) string =[]
              wanted_export (1,:) string =[]
           end
             obj.medium.manual_set_boundaries.unwanted_import = unwanted_import;
             obj.medium.manual_set_boundaries.unwanted_export = unwanted_export;
             obj.medium.manual_set_boundaries.wanted_import = wanted_import;
             obj.medium.manual_set_boundaries.wanted_export = wanted_export;
           
        end
    end

    end


    
function consistency_check(model)
    %UNTITLED4 Construct an instance of this class
    %   Detailed explanation goes here    

    disp("Check consistency of input model!")
    if isempty(model.rxns(find(contains(model.rxns,'biomass'))))
        error("You lost your objective function, when running fastcc!")
    end
    % check if the created model is now really consistent
    A = fastcc(model, 1e-4, 1);
    if length(model.rxns) ~= length(A)
        error("Your intput model is not consistent after running fastcc! Check your model!")
    else
        disp("Your model is consistent!")
    end
end

function mustBeMember_inModel(models_to_compare,obj)
            assert(all(ismember(models_to_compare,obj.run_names)), "The specified models to compare are not part of the model! Check the obj.run_names slot to get the correct labels!")
end






function [sampling_ordered] = get_sampling_orig_order(m,s,rxns_orig)
                                [~,mapping_rxns_in_orig_idx] = ismember(m.rxns,rxns_orig);
                                sampling_values = zeros(length(rxns_orig),size(s,2));
                                sampling_values(mapping_rxns_in_orig_idx,:) = s;
                                sampling_ordered = sampling_values;
%                                 id_biomass_ordered = find(matches(rxns_orig,"biomass_reaction"));
%                                 id_biomass = find(matches(m.rxns,"biomass_reaction"));
                                
%                                 min(sampling_values(id_biomass_ordered,:))
%                                 max(sampling_values(id_biomass_ordered,:))
%                                 
%                                 min(s(id_biomass,:))
%                                 max(s(id_biomass,:))
                                
end

function [sampling_fluxsum_ordered] = get_sampling_orig_order_mets(m,s,mets_all)
                                        [~,mapping_mets_in_orig_idx] = ismember(m.mets,mets_all);
                                        sampling_fluxsum_values = zeros(length(mets_all),size(s,2));
                                        sampling_fluxsum_values(mapping_mets_in_orig_idx,:) = s;
                                        sampling_fluxsum_ordered = sampling_fluxsum_values;
                               
end


function computeAndStoreFVAsim(i, rowIdx, colIdx, fvaData, exp)
    y = rowIdx(i);
    x = colIdx(i);

    [overallSim, rxnSim] = FVAsimilarity(fvaData{y}, fvaData{x});
    
    % Store symmetric results
    exp.fva_similarity(y, x) = overallSim;
    exp.fva_similarity(x, y) = overallSim;

    exp.fva_similarity_rxns{y, x} = rxnSim;
    exp.fva_similarity_rxns{x, y} = rxnSim;
end




function mustBeFileType(file_path,needed_file_format_ending)
    arguments
        file_path (1,1) string
        needed_file_format_ending (1,1) string
    end
            assert(exist(file_path,'file') ==2 , "Does the file exist ? Check again!")
            assert(~isempty(regexp(file_path,needed_file_format_ending + "$")),...
                   "Input must be a " + needed_file_format_ending + " file!!")
end





function scr_para = read_in_run_def_file(def_run_file)
    % Read the file
    C = readcell(def_run_file);
    
    % Separate field names and values (skip header)
    fields = C(2:end, 1);   % first column: slot_name
    values = C(2:end, 2);   % second column: value
    
    % Optional: convert numeric values stored as cell doubles to scalars
    for i = 1:numel(values)
        if iscell(values{i}) && isnumeric(values{i}) && numel(values{i})==1
            values{i} = values{i};  % already numeric, keep as-is
        elseif ismissing(values{i})
            values{i} = [];  % replace missing with empty
        elseif isstring(values{i}) || ischar(values{i})
            values{i} = string(values{i});  % ensure string
        end
    end
    
    % Create struct
    scr_para = cell2struct(values, fields, 1);
end



