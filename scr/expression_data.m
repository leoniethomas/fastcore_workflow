classdef expression_data
    %expressiond data -> performs the dscretization and enables visualization of the input sequencing data
    %   Detailed explanation goes here
    
    properties
        sample_names 
        feature_names_raw
        feature_names_norm
        raw_counts
        metadata
        normalized_counts
        discretized
    end
    
    methods
        function obj = expression_data(path_raw_counts,path_metadata,sample_label_column)
            if nargin ==3
                % read in sample names 
                obj.metadata = readtable(path_metadata, 'Delimiter','\t');
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
                
            end
        end
        function obj = get_normalized_data(obj, file_path)
            if nargin ==2
                obj.normalized_counts = readcell(file_path); 
            else
                error("fpkm normalization to be added to this function!")
            end
            
            % get sample and feature names 
            % bring the raw data into the right format
            find_numeric_entries = cellfun(@(x) isnumeric(x), obj.normalized_counts);
            sample_row = find(sum(find_numeric_entries,2)==0);
            feature_column = find(sum(find_numeric_entries)==0);
            obj.feature_names_norm = rmmissing(string(obj.normalized_counts(:,feature_column)));
            sample_names = rmmissing(string(obj.normalized_counts(sample_row,:)));
            obj.feature_names_norm(find(matches(obj.feature_names_raw,sample_names))) = [];
            obj.normalized_counts(sample_row,:) = [];
            obj.normalized_counts(:,feature_column) = [];
            % we extracted the feature and sample column
            % now check if the rest of the matrix is all numeric
            % then transfert to matrix
            if ~sum(sum(cellfun(@(x) ~isnumeric(x), obj.normalized_counts)) ~= 0)
                obj.normalized_counts =cell2mat(obj.normalized_counts);
            else
            end
            
            
            % check that the samples are in the corresponding order in
            % the metadata as well as in the expresion data           
            [a,b,c] = intersect(obj.sample_names, sample_names);
            if length(obj.sample_names) ~= length(a) | sum(c ~= b) > 0
                % if the number of samples is not equal or if the samples are not in the same order                    
                obj.sample_names = obj.sample_names(b);
                sample_names = sample_names(c);
                obj.metadata= obj.metadata(b,:);
                obj.raw_counts = obj.raw_counts(:,b);
                obj.normalized_counts = obj.normalized_counts(:,c);
            end

        end
        
        function obj = get_discretized_data(obj,figflag,file_path)
            mkdir(file_path + "Discretization" )
            obj.discretized = discretize_FPKM(obj.normalized_counts, obj.sample_names,figflag,char(file_path + "Discretization" + filesep));
            
            num_disc = hist(obj.discretized,3);
            figure
            %bar(perc_disc','stacked')
            bar(num_disc','stacked')
            xticks(1:length(obj.sample_names))
            xticklabels(obj.sample_names)
            xtickangle(90); 
            hold off
            saveas(gcf, file_path +  "discretized_count.png");
            
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
        
        function [coeff,score,latent,tsquared,explained,cluster] = QC_pca_kmeans(obj,data_slot,group_column, vis_pcs,num_k,vis_genes_idx,save_fig, pat_rep_label)
            
            data = obj.(data_slot);
            sample_name = string(obj.metadata.(group_column));
            
            [coeff,score,latent,tsquared,explained] = pca(data(vis_genes_idx,:)');
            cluster =num2cell(num2str(kmeans(data(vis_genes_idx,:)',num_k)));

            figure
            hold on
            for x = unique(sample_name)'
                %disp(x);
                idx = contains(sample_name,x{:});
                scatter(score(idx,vis_pcs(2)),score(idx,vis_pcs(1)))
            end   

            title('PCA - label')
            xlabel(['PC ', num2str(vis_pcs(2)), '  : ', num2str(explained(vis_pcs(2)))])
            ylabel(['PC ', num2str(vis_pcs(1)), '  : ', num2str(explained(vis_pcs(1)))])
            legend(regexprep(unique(sample_name)', pat_rep_label(1), pat_rep_label(2)) , ...
                   'location','best')
            hold off
            saveas(gcf, save_fig);

            figure
            hold on
            for x = unique(cluster)'
                %disp(x);
                idx = contains(cluster,x{:});
                scatter(score(idx,vis_pcs(2)),score(idx,vis_pcs(1)))
            end   

            title('PCA - cluster id')
            xlabel(['PC ', num2str(vis_pcs(2)), '  : ', num2str(explained(vis_pcs(2)))])
            ylabel(['PC ', num2str(vis_pcs(1)), '  : ', num2str(explained(vis_pcs(1)))])
            legend(unique(cluster), ...
                   'location','best')
            hold off
            saveas(gcf, save_fig);
        end
    end
end

























