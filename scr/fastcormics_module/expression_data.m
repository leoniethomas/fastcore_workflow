classdef expression_data
    % expression_data is a class defined to perform all the QC and
    % preprocessing steps needed to perform an analysis with fastcormics
    % The expression_data class provides a set of methods 
    %    -> expression_data initializting function -> reads in all the
    %       needed data 
    %       
    %
    
    properties
        sample_names % names of all samples/cells in the gene expression file read in
        metadata % metadata defining the properties of the different samples/cells
        feature_names_raw % defining the gene names used in the raw count data stored in the raw_counts slot
        raw_counts % raw unnormalized counts
        feature_names_norm % defining the gene names used for the normalized data slots 
        FPKM % Fragments Per Kilobase Million data
        TPM %  Transcripts Per Kilobase Million data
        vst_normalized_counts % variance stabilized data from Deseq2
        discretized % discretized data - generated useing fastcormics discretize_FPKM function
        features_metabolic_genes % features in the feature_names_norm slot which can be found in a metabolic model
        rxn_names % rxn names from the metabolic model
        mapping_exp_2_rxns % activity score for each rxn in rxn names per sample - generated using the fastcormics function map_expression_2_data_rFASTCORMICS
        pca
        source
    end
    
    methods
        function obj = expression_data(path_raw_counts,path_metadata,sample_label_column)
                % this function reads in the needed expression data to
                % perform the preprocessing & QC for fastcormics.
                %   path_raw_counts:     full path to the file storing the
                %                        counts
                %   path_metadata:       full path to the metadata file, defining
                %                        the characteristics of the samples/cells
                %                        in the expression data 
                %   sample_label_column: give the column after which the
                %                        samples are defined, this is
                %                        used later when visualizing the
                %                        data
                arguments
                   path_raw_counts (1,1) string {mustBeFileType(path_raw_counts,"txt")}
                   path_metadata (1,1) string {mustBeFileType(path_metadata,"txt")}
                   sample_label_column (1,1) string
                end
            
                
                % read in sample names 
                obj.metadata = readtable(path_metadata, 'Delimiter','\t');
                
                if ~any(contains(obj.metadata.Properties.VariableNames,sample_label_column))
                    disp("These are the columnnames in the metadata file:")
                    obj.metadata.Properties.VariableNames
                    error("The sample label column does not exist in the metadata! Check please which column you want to use to label the samples!")
                end
                obj.sample_names = string(obj.metadata.(sample_label_column))';
                obj.raw_counts = readcell(path_raw_counts);
                
                % bring the raw data into the right format
                find_numeric_entries = cellfun(@(x) isnumeric(x), obj.raw_counts);
                sample_row = find(sum(find_numeric_entries,2)==0);
                feature_column = find(sum(find_numeric_entries)==0);
                obj.feature_names_raw = rmmissing(string(obj.raw_counts(:,feature_column)));
                sample_names = rmmissing(string(obj.raw_counts(sample_row,:)));
                obj.feature_names_raw(find(matches(obj.feature_names_raw,sample_names))) = [];
                obj.raw_counts(sample_row,:) = [];
                obj.raw_counts(:,feature_column) = [];
                % we extracted the feature and sample column
                % now check if the rest of the matrix is all numeric
                % then transfert to matrix
                if ~sum(sum(cellfun(@(x) ~isnumeric(x), obj.raw_counts)) ~= 0)
                    obj.raw_counts =cell2mat(obj.raw_counts);
                else
                    error("The raw count matrix entails other non-numeric columns/rows than the features and the samples! Check your matrix before using it as input!")
                end
              
                if size(obj.raw_counts,2) ~= length(sample_names)
                    error("The choosen sample name row does not entail all the sample ids, to cover the number of columns in the raw count data. Check that there are as many identifiers in the sample row as there are columns in the data!")
                end
                
                % check that the samples are in the corresponding order in
                % the metadata as well as in the expresion data           
                [a,b,c] = intersect(obj.sample_names, sample_names);
                if length(obj.sample_names) ~= length(a) | sum(c ~= b) > 0
                    % if the number of samples is not equal or if the samples are not in the same order                    
                    obj.sample_names = obj.sample_names(b);
                    sample_names = obj.samples(c);
                    obj.metadata= obj.metadata(b,:);
                    obj.raw_counts = obj.raw_counts(:,c);
                end
                disp("expression data object created:")
                obj
                disp("expression data entails " + string(length(obj.sample_names)) + " samples!")
            
        end
        function obj = get_normalized_data(obj, file_path, slot)
            % This function reads in the normalized data 
            % normalized data referres here to TPM and/or FPKM
            % (vst_normalized counts)
            
            arguments
               obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
               file_path (1,1) string {mustBeFileType(file_path,"csv")}
               slot (1,1) string {mustBeMember(slot,["TPM","FPKM","vst_normalized"])}
            end

            obj.(slot) = readcell(file_path);          
            % get sample and feature names 
            % bring the raw data into the right format
            find_numeric_entries = cellfun(@(x) isnumeric(x), obj.(slot));
            sample_row = find(sum(find_numeric_entries,2)==0);
            feature_column = find(sum(find_numeric_entries)==0);
            features_in_data = rmmissing(string(obj.(slot)(:,feature_column)));
            sample_names = rmmissing(string(obj.(slot)(sample_row,:)));
            obj.(slot)(sample_row,:) = [];
            obj.(slot)(:,feature_column) = [];
            
            features_in_data(find(matches(features_in_data,sample_names))) = [];
            if ~isempty(obj.feature_names_norm)
                % in case there has been normalized data already read in
                % the new data is filtered according what is already stored
                % in the feature_name_norm 
                [~,idx_in_feature_names_norm] = ismember(obj.feature_names_norm, features_in_data);
                obj.(slot) = obj.(slot)(idx_in_feature_names_norm,:);
                if length(idx_in_feature_names_norm) < length(features_in_data)
                   disp("You lost some of your features in the " + slot + " data slot since it was indexed based on the features in the .feature_names_norm slot of the expression_data object!") 
                end
            else
                obj.feature_names_norm = features_in_data;
            end
            
            
            
            % we extracted the feature and sample column
            % now check if the rest of the matrix is all numeric
            % then transfer to matrix
            if ~sum(sum(cellfun(@(x) ~isnumeric(x), obj.(slot))) ~= 0)
                obj.(slot) =cell2mat(obj.(slot));
            else
            end
            
            
            % check that the samples are in the corresponding order in
            % the metadata as well as in the expresion data   
            if length(sample_names) ~= size(obj.(slot),2)
                error("There are more entries in the sample row than there are data in the data, check which of the entries you need to get rid of.")
            end
            [~,b] = intersect(obj.sample_names, sample_names);
            obj.(slot) = obj.(slot)(:,b);
            if length(b) ~= length(obj.sample_names)
               error("Not all the samples from the metadata could be found in the data!") 
            end
            
            disp("expression data object created:")
            obj
            disp("expression data entails " + string(length(obj.sample_names)) + " samples!")

        end
        
        function obj = map_expression_2_rxns(obj,model_used,dic_gene_ids_entrez_used)
           % this function executes the mapping of the gene expression to the rnx in the model
            arguments
               obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
               model_used (1,1) struct
               dic_gene_ids_entrez_used (:,:) table
            end
            
            if exist("map_expression_2_data_rFASTCORMICS",'file') == 0
               error("Fastcormics is not installed or the installation was not added to the path variable! The function used for mapping can not be found! map_expression_2_data_rFASTCORMICS function is needed to execute this task!") 
            end
           
            mapping = map_expression_2_data_rFASTCORMICS(model_used, ...
                                                         obj.discretized,...
                                                         dic_gene_ids_entrez_used,...
                                                         obj.feature_names_norm);
            mapping = sparse(mapping);
            obj.mapping_exp_2_rxns = mapping;
            obj.rxn_names = string(model_used.rxns);
        end
        
        function obj = get_discretized_data(obj,figflag,file_path_results,slot)
            % This function executes the discretize function of fastcormics
            % and saves the output to a folder
            arguments
               obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
               figflag (1,1) double {mustBeMember(figflag,[1,0])}
               file_path_results (1,1) string
               slot (1,1) string {mustBeMember(slot,["TPM","FPKM","vst_normalized","raw_counts"])} ="FPKM"
            end

            mkdir(file_path_results + "Discretization" )
            obj.discretized = discretize_FPKM(obj.(slot), ...
                                              obj.sample_names,figflag,...
                                              char(file_path_results + "Discretization" + filesep));
            
            num_disc = hist(obj.discretized,3);
            figure
            %bar(perc_disc','stacked')
            bar(num_disc','stacked')
            xticks(1:length(obj.sample_names))
            xticklabels(obj.sample_names)
            xtickangle(90); 
            hold off
            saveas(gcf,file_path_results +  "discretized_count.png");
            
        end
        
        function obj = get_metabolic_genes(obj,model_used,dic_gene_ids_entrez_used)
            % 
            
            arguments
               obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
               model_used (1,1) struct
               dic_gene_ids_entrez_used table
            end        
            
            metabolic_genes_entrez = string(regexprep(model_used.model.genes,".1",""));
            obj.features_metabolic_genes = string(dic_gene_ids_entrez_used.Var1.dico.SYMBOL(find(matches(dic_gene_ids_entrez_used.Var1.dico.ENTREZ,...
                                                                   metabolic_genes_entrez))));
            
        end
        
        function obj = get_QC_plots(obj,data_type,group_column,file_path)
            plot_data = obj.(data_type);
            sample_names = string(obj.metadata.(group_column));
            
            bar(sum(plot_data))
            title("Number of reads per sample with " + data_type + " : ")
            xlabel("samples")
            ylabel("# of reads")
            xticks(1:length(sample_names))
            xticklabels(sample_names)
            xtickangle(90); 
            hold off
            saveas(gcf,file_path +  data_type + "_barplot.png");

            bar(sum(plot_data == 0,1))
            title("Number of 0s per sample with " + data_type + " : ")
            xlabel("samples")
            ylabel("# of zeros")
            xticks(1:length(sample_names))
            xticklabels(sample_names)
            xtickangle(90); 
            hold off
            saveas(gcf,file_path + data_type +"_zero_per_sample_barplot.png");
            
            figure
            boxplot(plot_data)
            title("expression per sample -> with " + data_type + " data")
            xlabel("samples")
            ylabel("expression in counts")
            xticks(1:length(sample_names))
            xticklabels(sample_names)
            xtickangle(90); 
            hold off
            saveas(gcf, file_path + data_type + "_boxplots.png");
            
            data2=plot_data;
            data2(data2==0)=NaN;
            figure
            boxplot(data2)
            title("expression per sample without 0 -> with " + data_type + " data")
            xlabel("samples")
            ylabel("expression in counts")
            xticks(1:length(sample_names))
            xticklabels(sample_names)
            xtickangle(90); 
            hold off
            saveas(gcf, file_path + data_type + "_boxplots_nozeros.png");
            
            figure
            hold on
            fpkm = plot_data;
            fpkm(fpkm==0)=NaN; %remove zeros for densityplot
            fpkm=log2(fpkm); %log2 scaling
            for i=1:size(fpkm,2)
                [probability_estimate,xi] = ksdensity(fpkm(:,i));
                plot(xi,probability_estimate,':k','LineWidth',1);
            end
            title("Distribution per sample with " + data_type + " data")
            hold off
            saveas(gcf, file_path + data_type + "_density_dist.png");

            
            
        end
        
        function obj = perform_pca_kmeans(obj,data_slot,num_k,vis_features)
            
            arguments
               obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
               data_slot (1,1)
               num_k (1,1)
               vis_features (1,:) =ones(1,size(obj.(data_slot),1))
            end  
            
            
            data = obj.(data_slot);
            data = full(data(find(vis_features),:));
            
            [obj.pca.(data_slot).coeff,obj.pca.(data_slot).score,obj.pca.(data_slot).latent,obj.pca.(data_slot).tsquared,obj.pca.(data_slot).explained] = pca(data');
            cluster =num2cell(num2str(kmeans(data',num_k)));
            column_name = "kmeans_k" + string(num_k) + "_" + data_slot + "_features_" + string(length(find(vis_features)));
            obj.metadata.(column_name) = str2num(cell2mat(cluster));
            
        end
        
        function visualize_dimreduction(obj,colour_label,reduction,data_slot,save_fig, pat_rep_label,vis_dim)

            arguments
               obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
               colour_label (1,1) ="Treatment"
               reduction (1,1) ="pca"
               data_slot (1,1) ="FPKM"
               save_fig (1,1) ="mapped_expression_" + reduction + ".png"
               pat_rep_label (1,2) =["_" "/_"]
               vis_dim (1,2) =[1,2]
            end  
            
            sample_cat = string(obj.metadata.(colour_label));
            score = obj.(reduction).(data_slot).score;
            explained = obj.(reduction).(data_slot).explained;
            figure
            hold on
            for x = unique(sample_cat)'
                %disp(x);
                idx = contains(sample_cat,x{:});
                scatter(score(idx,vis_dim(2)),score(idx,vis_dim(1)))
            end   

            title(regexprep(reduction+" - " + data_slot + " - " + colour_label, "_" , " "))
            xlabel(['PC ', num2str(vis_dim(2)), '  : ', num2str(explained(vis_dim(2)))])
            ylabel(['PC ', num2str(vis_dim(1)), '  : ', num2str(explained(vis_dim(1)))])
            legend(regexprep(unique(sample_cat)', pat_rep_label(1), pat_rep_label(2)) , ...
                   'location','best')
            hold off
            saveas(gcf, regexprep(save_fig,"PCA.png","PCA_label.png"));
        end
    end
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


function mustBeValid_expression_data_object(object)
    
    % check the number of samples in the different slots of the expression
    % object!
    assert(length(object.sample_names) == size(object.metadata,1), ...
           "Your sample_names slot and the metadata do not hold the same number of samples! Check your object before supplying it to the function!")
    
    data_slots_to_be_checked = ["FPKM" "TPM" "vst_normalized_counts" "raw_counts" "discretized" "mapping_exp_2_rxns"];
    data_slot_size = arrayfun(@(x) size(object.(x),2),...
                              data_slots_to_be_checked);
    correct_size = data_slot_size == 0 | data_slot_size == length(object.sample_names);
    
    assert(all(correct_size),...
           "At least one of your slots of the expression_data does not have as much samples as there are listed in the samples slot! Check these slots: " + strjoin(data_slots_to_be_checked(find(~correct_size)), ", "));
    
    % check that the dimension of rxn names and mapping data matrix aggree
    assert(size(object.mapping_exp_2_rxns,1) == length(object.rxn_names) | isempty(object.mapping_exp_2_rxns), ...
          "The rxns names specified in the rxn_names slot does not correspond to the number of rows in the mapping_exp_2_rxns slot! Check your object!!")
    
    % check that the raw counts rows aggree with the dimension the raw data matrix has!
    assert(size(object.raw_counts,1) == length(object.feature_names_raw) ,...
           "The gene names saved in the fature_names_raw slot do not have the same length compared to the matrix that is stored in the raw_counts slot! Check! These two must aggree!")
    
    % check the aggreement of the features stored in the
    % features_names_norm slot 
    data_norm_features_count = ["FPKM" "TPM" "discretized" "vst_normalized_counts"];
    data_slot_size = arrayfun(@(x) size(object.(x),1),...
                              data_norm_features_count);
    correct_feature_count = data_slot_size == 0 | data_slot_size == length(object.feature_names_norm);
    assert(all(correct_feature_count),...
           "At least one of your slots of the expression_data does not have as much genes as there are listed in the feature_names_norm slot! Check these slots: " + strjoin(data_norm_features_count(find(~correct_feature_count)), ", "));
    
end


    











