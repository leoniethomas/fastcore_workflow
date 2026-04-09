classdef expression_data
    % expression_data is a class defined to perform all the QC and
    % preprocessing steps needed to perform an analysis with fastcormics
    % The expression_data class provides a set of methods 
    %    -> expression_data initializting function -> reads in all the
    %       needed data 
    %       
    %
    
    properties
        Properties
        sample_names % names of all samples/cells in the gene expression file read in
        metadata % metadata defining the properties of the different samples/cells
        raw_counts % raw unnormalized counts
        norm_counts %counts normalized for library size/gene length(for non umi RNAseq)
        znorm_counts % scaled and centered FPKM
        discretized % discretized data - generated useing fastcormics discretize_FPKM function
        features_metabolic_genes % features in the feature_names_norm slot which can be found in a metabolic model
        mapping_exp_2_rxns % activity score for each rxn in rxn names per sample - generated using the fastcormics function map_expression_2_data_rFASTCORMICS
        model_id
        FPKM
        feature_names_norm
        dico
        pca
        umap
    end
    
    methods
        function obj = expression_data(path_raw_counts,path_metadata,dico,sample_label_column)
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
                % arguments
                %    path_raw_counts (1,1) %string {mustBeFileType(path_raw_counts,"csv")}
                %    path_metadata (1,1) %string {mustBeFileType(path_metadata,"txt")}
                %    sample_label_column (1,1) string
                % end
            
                
                % read in sample names 
                obj.metadata = readtable(path_metadata);
                obj.dico = dico;
                
                if ~any(contains(obj.metadata.Properties.VariableNames,sample_label_column))
                    disp("These are the columnnames in the metadata file:")
                    obj.metadata.Properties.VariableNames
                    error("The sample label column does not exist in the metadata! Check please which column you want to use to label the samples!")
                end
                obj.Properties = string(properties(obj));
                obj.sample_names = string(obj.metadata.(sample_label_column))';
                obj.raw_counts = readtable(path_raw_counts);
                
                % bring the raw data into the right format
                feature_column = find(obj.raw_counts.Properties.VariableTypes ~= "double");   
                % add postfix to the genes which are double in the table 
                gene_names = string(obj.raw_counts{:,feature_column});
                if length(gene_names) ~= length(unique(gene_names))
                    gene_names = addpostfixtogeneidentifier(gene_names);
                end
                if isempty(obj.raw_counts.Properties.RowNames)
                   obj.raw_counts.Properties.RowNames = gene_names;
                   obj.raw_counts(:,feature_column) = [];
                end
                if isempty(obj.raw_counts.Properties.VariableNames)
                    error("The Variablenames are empty for the raw counts table, adjust code to get the sample names as Column names!!")
                end
              
                % sort the table according to the metadata given 
                if sum(matches(obj.sample_names, string(obj.raw_counts.Properties.VariableNames))) == 0
                    obj.sample_names = regexprep(obj.sample_names,"_", "-");
                    obj.raw_counts.Properties.VariableNames = regexprep(obj.raw_counts.Properties.VariableNames,"_", "-");
                    if sum(~matches(obj.sample_names, string(obj.raw_counts.Properties.VariableNames))) ~= 0
                        error("Since non of the sample names between the metadata matched _ was replaced by - but still there are samples in the metadata that can not be found in the raw count data! Check your csv files!")
                    end
                end

                obj.raw_counts = obj.raw_counts(:,obj.sample_names);
            
        end
        function obj = get_normalized_data(obj, file_path)
            % This function reads in the normalized data 
            % normalized data referres here to TPM and/or FPKM
            % (vst_normalized counts)
            
            % arguments
            %    obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
            %    file_path (1,1) string {mustBeFileType(file_path,"csv")}
            % end

            obj.norm_counts = readtable(file_path);

            % bring the raw data into the right format
            feature_column = find(obj.norm_counts.Properties.VariableTypes ~= "double"); 

            % add postfix to the genes which are double in the table 
            gene_names = string(obj.norm_counts{:,feature_column});
            if length(gene_names) ~= length(unique(gene_names))
               gene_names = addpostfixtogeneidentifier(gene_names);
            end
            if isempty(obj.norm_counts.Properties.RowNames)
               obj.norm_counts.Properties.RowNames = gene_names;
               obj.norm_counts(:,feature_column) = [];
            end
            if isempty(obj.norm_counts.Properties.VariableNames)
                error("The Variablenames are empty for the normalized counts table, adjust code to get the sample names as Column names!!")
            end
          
            % sort the table according to the metadata given 
            if sum(matches(obj.sample_names, string(obj.norm_counts.Properties.VariableNames))) == 0
                obj.sample_names = regexprep(obj.sample_names,"_", "-");
                obj.norm_counts.Properties.VariableNames = regexprep(obj.norm_counts.Properties.VariableNames,"_", "-");
                if sum(~matches(obj.sample_names, string(obj.norm_counts.Properties.VariableNames))) ~= 0
                    error("Since non of the sample names between the metadata matched _ was replaced by - but still there are samples in the metadata that can not be found in the raw count data! Check your csv files!")
                end
            end

            obj.norm_counts = obj.norm_counts(:,obj.sample_names);

        end
        
        function obj = map_expression_2_rxns(obj,model_used,dic_gene_ids_entrez_used)
           % this function executes the mapping of the gene expression to the rnx in the model
            % arguments
            %    obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
            %    model_used (1,1) struct
            %    dic_gene_ids_entrez_used (:,:) table
            % end
            
            if exist("mapExpressionToModel",'file') == 0
               error("Fastcormics is not installed or the installation was not added to the path variable! The function used for mapping can not be found! map_expression_2_data_rFASTCORMICS function is needed to execute this task!") 
            end
           
            mapping = mapExpressionToModel(model_used, ...
                                                         obj.discretized,...
                                                         dic_gene_ids_entrez_used,...
                                                         obj.norm_counts.Properties.RowNames,1);
            obj.model_id = model_used.id;
            mapping = sparse(mapping);
            obj.mapping_exp_2_rxns = array2table(mapping, ...
                                             'RowNames', model_used.rxns, ...
                                             'VariableNames', obj.norm_counts.Properties.VariableNames);
        end
        
        function obj = get_discretized_data(obj,figflag,file_path_results,fig_format)
            % This function executes the discretize function of fastcormics
            % and saves the output to a folder
            arguments
               obj 
               figflag (1,1) 
               file_path_results (1,1) string
               fig_format ='.png'
            end

            mkdir(file_path_results)
            [obj.discretized, ...
             obj.znorm_counts] = discretizeFPKM(obj.norm_counts{:,:}, ...
                                                   obj.sample_names,figflag,...
                                                   char(file_path_results), fig_format);
            
            vals = obj.discretized(:);              % flatten if needed
            [uVals, ~, idx] = unique(vals);        % unique values
            
            % Count occurrences per sample (column-wise)
            num_disc = zeros(numel(uVals), size(obj.discretized,2));
            for i = 1:numel(uVals)
                num_disc(i,:) = sum(obj.discretized == uVals(i), 1);
            end
            
            figure
            bar(num_disc', 'stacked')
            
            xticks(1:length(obj.sample_names))
            xticklabels(obj.sample_names)
            xtickangle(90)
            
            % Legend automatically from unique values
            legend(string(uVals), 'Location','best')
            hold off
            saveas(gcf,char(file_path_results ) +  "discretized_count.svg");
            saveas(gcf,char(file_path_results ) +  "discretized_count.png");

        end
        
        function obj = get_metabolic_genes(obj,model_used, wanted_id)
            % 
            
            % arguments
            %    obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
            %    model_used (1,1) struct
            %    dic_gene_ids_entrez_used table
            % end        
            
            metabolic_genes_entrez = regexprep(string(model_used.(string(fieldnames(model_used))).genes), '\.1$', '');    
            % find out which type of gene id is used by scanning the dico for matches 
            f = @(x) sum(matches(metabolic_genes_entrez,string(obj.dico{:,x})));
            match_per_column_dico = arrayfun(f, 1:size(obj.dico,2));
            [~,mapping_column_idx] = max(match_per_column_dico);
            mapping_column = obj.dico{:,mapping_column_idx};

            find_matches_in_data = find(matches(mapping_column,metabolic_genes_entrez));
            
            obj.features_metabolic_genes = string(obj.dico.(wanted_id));
            obj.features_metabolic_genes = obj.features_metabolic_genes(find_matches_in_data);
            
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
        
        function obj = delete_sample(obj,sample_names)
            arguments
                obj (1,1) expression_data  {mustBeValid_expression_data_object(obj)}
                sample_names (1,:) string
            end
            
            sample_count_init = length(obj.sample_names);
            keep_sample_idx = find(~ismember(obj.sample_names,sample_names));
            obj.sample_names = obj.sample_names(keep_sample_idx);
            obj.metadata = obj.metadata(keep_sample_idx,:);

            slots_with_sample_data = obj.Properties(find(sample_count_init == cellfun(@(p) size(obj.(p),2), obj.Properties)));
            for slot = slots_with_sample_data'
                obj.(slot) = obj.(slot)(:,keep_sample_idx);
            end

            % TODO -> delete all the dimension reductions already
            % performed, need to be performed from scratch, since samples
            % were removed
            
        end
        function obj = perform_pca_kmeans(obj,data_slot,num_k,perform_clustering_on_pca,vis_features)
            
            arguments
               obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
               data_slot (1,1)
               num_k (1,1)
               perform_clustering_on_pca (1,1) =1
               vis_features (1,:) =ones(1,size(obj.(data_slot),1))
            end  
            
            
            data = obj.(data_slot);
            data = full(data(find(vis_features),:));
            
            
            [obj.pca.(data_slot).coeff,obj.pca.(data_slot).score,obj.pca.(data_slot).latent,obj.pca.(data_slot).tsquared,obj.pca.(data_slot).explained] = pca(data');
            if perform_clustering_on_pca
                cluster =num2cell(num2str(kmeans(obj.pca.(data_slot).score,num_k)));
            else
                cluster =num2cell(num2str(kmeans(data',num_k)));
            end
            column_name = "kmeans_k" + string(num_k) + "_" + data_slot + "_features_" + string(length(find(vis_features)));
            obj.metadata.(column_name) = str2num(cell2mat(cluster));
            
        end
        
        function obj = perform_umap(obj,data_slot,num_neighbors,vis_features)
            
            arguments
               obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
               data_slot (1,1)
               num_neighbors (1,1) =3
               vis_features (1,:) =ones(1,size(obj.pca.(data_slot).score,1))
            end  
            
            
            data_plot = obj.pca.(data_slot).score;
            data_plot = full(data_plot(:,find(vis_features)));
            % set seed 
            rng(123);
            % --- Run UMAP dimensionality reduction ---
            % Make sure run_umap.m is on your MATLAB path
            [reduction, ~] = run_umap( ...
                data_plot, ...
                'n_components', 2, ...     % project to 2D
                'n_neighbors', num_neighbors, ...     % size of local neighborhood
                'min_dist', 0.1, ...       % how tightly points are packed
                'metric', 'euclidean', ...  % distance metric
                'verbose', false, ...
                'random_state', 42 ... % makes sure that the plot is reproducible
            );
            obj.umap.(data_slot).score = reduction;
            
        end

        function save_data_as_csv(obj,slot,path,rowname_slot,colname_slot)
            arguments
                obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
                slot (1,1) string ="FPKM"
                path (1,1) string = "./saved_data.csv"
                rowname_slot (1,1) string ="feature_names_norm"
                colname_slot (1,1) string ="sample_names"
            end
            
            T = array2table(obj.(slot), 'VariableNames', obj.(colname_slot), 'RowNames', obj.(rowname_slot));
            writetable(T, path, 'WriteRowNames', true);
        end

        function obj = fix_duplicated_feature_names(obj,slot)

            % Original row names
            rownames = obj.(slot);
            
            % Make row names unique by appending numeric suffixes
            [uniqueNames, ~, ic] = unique(rownames, 'stable'); 
            counts = accumarray(ic, 1);
            suffixes = zeros(size(rownames));
            
            for i = 1:length(rownames)
                if counts(ic(i)) > 1
                    suffixes(i) = sum(strcmp(rownames(1:i), rownames{i}));
                end
            end
            
            % Add postfix only to duplicates
            for i = 1:length(rownames)
                if suffixes(i) > 1
                    rownames{i} = [rownames{i} '_' num2str(suffixes(i))];
                end
            end
            obj.(slot) = rownames;

        end
        
        function visualize_dimreduction(obj,colour_label,reduction,data_slot,shape_label,vis_dim, save_fig, pat_rep_label)

            arguments
               obj (1,1) expression_data {mustBeValid_expression_data_object(obj)}
               colour_label (1,1) ="Treatment"
               reduction (1,1) ="pca"
               data_slot (1,1) ="FPKM"
               shape_label (1,1) = "None"
               vis_dim (1,2) =[1,2]
               save_fig (1,1) ="mapped_expression_" + reduction + ".png"
               pat_rep_label (1,2) =["_" " "]
            end  
            
            sample_cat = string(obj.metadata.(colour_label));
            if shape_label == "None"
                shape_cat  = repmat("o",length(sample_cat),1);
            else
                shape_cat  = string(obj.metadata.(shape_label));   % NEW
            end

            score = obj.(reduction).(data_slot).score;
            
            figure
            ax1 = axes; hold on
            
            % --- Plot scatter points ---
            cats_color = unique(sample_cat)';
            
            cats_shape = unique(shape_cat)';
            markers = {'o','s','^','d','v','p','h','x','+'};
            colors = lines(numel(cats_color));
            
            for i = 1:numel(cats_color)
                for j = 1:numel(cats_shape)
                    idx = sample_cat == cats_color(i) & shape_cat == cats_shape(j);
                    if any(idx)
                        scatter(score(idx,vis_dim(2)), score(idx,vis_dim(1)), ...
                                80, colors(i,:), 'filled', 'Marker', markers{mod(j-1,numel(markers))+1});
                    end
                end
            end
            
            title(regexprep(reduction+" - " + data_slot, pat_rep_label(1) , pat_rep_label(2)))
            
            % --- Color legend ---
            lgd1 = legend(arrayfun(@(i) plot(NaN,NaN,'o','MarkerFaceColor',colors(i,:),...
                                'MarkerEdgeColor',colors(i,:),'LineStyle','none','MarkerSize',8), ...
                                1:numel(cats_color)), cats_color, 'Location','northeast');
            ldg1.Color = "none";
            title(lgd1,regexprep(colour_label, pat_rep_label(1) , pat_rep_label(2)))
            
            % --- Shape legend (overlay axes) ---
            if length(unique(cats_shape)) ~= 1
            ax2 = copyobj(ax1,gcf); delete(ax2.Children); hold(ax2,'on')
            lgd2 = legend(arrayfun(@(j) plot(ax2,NaN,NaN,'Marker',markers{mod(j-1,numel(markers))+1},...
                                'LineStyle','none','MarkerSize',8,'Color','k'), 1:numel(cats_shape)), ...
                          cats_shape, 'Location','northwest');
            set(ax2,'Color','none','XTick',[],'YTick',[],'Box','off','Visible','off')
            title(lgd2,regexprep(shape_label, pat_rep_label(1) , pat_rep_label(2))); set(lgd2,'Color','none')

            set(lgd2,'Color','none')  % transparent background
            end

            
            if reduction == "pca"
                explained = obj.(reduction).(data_slot).explained;
                xlabel(['PC ', num2str(vis_dim(2)), '  : ', num2str(explained(vis_dim(2)))])
                ylabel(['PC ', num2str(vis_dim(1)), '  : ', num2str(explained(vis_dim(1)))])
            else
                xlabel([reduction,'', num2str(vis_dim(2))])
                ylabel([reduction,'', num2str(vis_dim(1))])
            end
      
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


    

function gene_names = addpostfixtogeneidentifier(gene_names)
    
    % Find groups of identical names
    [uniqueNames, ~, idx] = unique(gene_names, 'stable');
    
    % Count occurrences within each group
    counts = zeros(size(gene_names));
    for i = 1:numel(uniqueNames)
        mask = (idx == i);
        counts(mask) = 1:nnz(mask);
    end
    
    % Append suffix only to duplicates
    dupMask = counts > 1;
    gene_names(dupMask) = gene_names(dupMask) + "." + (counts(dupMask)-1);

end







