function project = visualize_sampling_landscape(project,comparison_name, rxn_to_visualize,options)
    % This function will visualize your reaction of interest on the a
    % dimension reduced space. 
    % 
    % Inputs: 
    %   - project:          the object which is the output of the single_model_analysis
    %                       entailing the results of fba,fva,sampling, single gene
    %                       deletion etc. for a single model 
    %   - comparison_name:  name of the comparsion the dimension reduction
    %                       should be performed
    %   - rxn_to_visualize: reaction values to be visualized as color on in the figure
    %   - options:          
    arguments
        project 
        comparison_name
        rxn_to_visualize (:,1) string =["biomass_reaction"]
        options.dim_reduction_type (1,1) string {mustBeMember(options.dim_reduction_type,["PCA", "UMAP"])} ="UMAP"
        options.pcs_vis (1,2) =[1,2]
        options.sampling_feature (1,1) string {mustBeMember(options.sampling_feature,["flux", "fluxsum"])} ="flux"
        options.num_clusters =0
        options.pcs_used_dim_red =0
        options.perform_kmeans =1
        options.thinning =5
        options.n_neighbors =40
        options.overwrite =0
    end
    if options.num_clusters == 0
        options.num_clusters = length(unique(project.comparisons.(comparison_name).sample_model_labels));
    end
    reference_model = project.comparisons.(comparison_name).reference_model;
    
    sample_model_labels = project.comparisons.(comparison_name).sample_model_labels;
    if options.sampling_feature == "flux"
        ordered_samples = project.comparisons.(comparison_name).ordered_samples;
        dim_names = project.models.(reference_model).model.rxns;
    elseif options.sampling_feature == "fluxsum"
        ordered_samples = project.comparisons.(comparison_name).ordered_samples_fluxsum;
        dim_names = project.models.(reference_model).model.mets;
    else
        error("You choose a  non valid argument for the sampling feature! Choose flux or fluxsum, by default the flux per reaction will be portrayed!")
    end
    
    % create new slot for output of function 
    project.comparisons.(comparison_name).dimension_reduction = struct();

    % run the Principle component analysis
    pc_x = options.pcs_vis(1);
    pc_y = options.pcs_vis(2);
    rxn_color = rxn_to_visualize;

    % filter out dimensions which are always 0
    ordered_samples_filter = ordered_samples(find(any(ordered_samples ~= 0,2)),:);
    ordered_dim_names = dim_names(find(any(ordered_samples ~= 0,2)),:);
    % z-scale the matrix so that all features have the same influence on
    % the dimension reduction
    samples_matrix = ordered_samples_filter';
    scaled_samples = zscore(samples_matrix, 0, 1); 
    % peform pca
    [pca_samp.coeff,pca_samp.score,pca_samp.latent,pca_samp.tsquared,pca_samp.explained] = pca(scaled_samples);
    % determine number of PCs used in downstream analysis or use the number
    % of PCs specified in the options parameter as function input
    if options.pcs_used_dim_red == 0
        cumulativeVariance = cumsum(pca_samp.explained);
        numPCs = find(cumulativeVariance >= 70, 1, 'first');
        fprintf('Number of PCs to reach 70%% variance: %d\n', numPCs);
    else
        numPCs = options.pcs_used_dim_red;
    end
    project.comparisons.(comparison_name).dimension_reduction.pca = pca_samp;
    
    
    % perform umap 
    if options.dim_reduction_type == "UMAP"
        % prepare data for dimension reduction 
        X = pca_samp.score(:,1:numPCs);
        thin = options.thinning;
        X = X(1:thin:end,:);
        fprintf('Computing umap dimension reduction ...!');
        [reduction,~] = run_umap(X, ...
                                'n_neighbors',options.n_neighbors, ...
                                'min_dist',0.3, ...
                                'verbose', 'none',...  
                                'method', 'C vectorized'); % 'MATLAB vectorized' also works, but much slower! C++ did not work for me
                                % 'MEX' is fastest but needs manual download of files for every execution
        fprintf('...finished!');
        %figure
        %scatter(reduction(:,1),reduction(:,2),20,categorical(labels),'filled')

        project.comparisons.(comparison_name).dimension_reduction.umap.reduction = reduction;
        project.comparisons.(comparison_name).dimension_reduction.umap.n_neighbors = options.n_neighbors;
        project.comparisons.(comparison_name).dimension_reduction.umap.thinning = thin;
        project.comparisons.(comparison_name).dimension_reduction.umap.used_samples_idx = 1:thin:length(sample_model_labels);
        project.comparisons.(comparison_name).dimension_reduction.umap.used_pcs = numPCs;
    end
    if options.perform_kmeans ==1
        % prepare data for dimension reduction 
        X = pca_samp.score(:,1:numPCs);
        thin = options.thinning;
        X = X(1:thin:end,:);
        labels = sample_model_labels(1:thin:end);
        kmeans_results = struct();
        [kmeans_results.idx, kmeans_results.C,...
         kmeans_results.sumd,kmeans_results.D] = kmeans(X, options.num_clusters, ...
                                                        'Distance','sqeuclidean', ...
                                                        'Replicates',10, ...
                                                        'MaxIter',500, ...
                                                        'Display','final');

        % figure;
        % scatter(reduction(:,1), reduction(:,2), 20, categorical(idx), 'filled');
        % xlabel('UMAP 1'); ylabel('UMAP 2');
        % title('K-means clusters on PCA space visualized via UMAP');
        % colorbar;

        % check quality of clustering - silhoutte score and homogeneity of
        % the clustering
        [kmeans_results.silhouette,~] = silhouette(X,kmeans_results.idx,'Euclidean', 'MaxIter', 100);
        mean_sil = mean(kmeans_results.silhouette);
        %
        homogen = [];
        for k = unique(kmeans_results.idx)'
            labels_in_cluster = labels(find(kmeans_results.idx == k));
            [s,~,j]=unique(labels_in_cluster);
            f = s{mode(j)};
            m = sum(labels_in_cluster == f);
            homogen = [ homogen, m/length(labels_in_cluster)];
        end
        kmeans_results.homogen = mean(homogen);

        project.comparisons.(comparison_name).dimension_reduction.kmeans = kmeans_results; 
    end
    
    % Visualization of computed metrics
    
    if options.dim_reduction_type == "PCA"

         pca_samp = project.comparisons.(comparison_name).dimension_reduction.pca;
         
         dim1 = pca_samp.score(:,pc_x);
         dim2 = pca_samp.score(:,pc_y);
         labels = sample_model_labels;
         
        
         figure
         scatter(dim1,dim2,20,...
                 categorical(labels),'filled')
         add_labels(dim1, dim2, labels)
        
         title("PC "+ num2str(pc_x) +  " & " +  num2str(pc_y) + " with the model lables!", 'FontSize',18)
         xlabel("PC" + num2str(pc_x) + " var: " + pca_samp.explained(pc_x), 'FontSize',18)
         ylabel("PC" + num2str(pc_y) + " var: " + pca_samp.explained(pc_y), 'FontSize',18)
         hold off
        

    else

        umap_res = project.comparisons.(comparison_name).dimension_reduction.umap;
         
        dim1 = umap_res.reduction(:,1);
        dim2 = umap_res.reduction(:,2);
        labels = sample_model_labels(umap_res.used_samples_idx);
        dim_used = umap_res.used_pcs;
         
        
        figure
        scatter(dim1,dim2,20,...
                 categorical(labels),'filled')
        add_labels(dim1, dim2, labels)
        
        title("UMAP with the model lables", 'FontSize',18)
        xlabel("UMAP1 " + num2str(dim_used) + "PCs -> 2D", 'FontSize',14)
        ylabel("UMAP2 " + num2str(dim_used) + "PCs -> 2D", 'FontSize',14)

        hold off

    end

    if options.perform_kmeans

        if options.dim_reduction_type == "PCA"
             pca_samp = project.comparisons.(comparison_name).dimension_reduction.pca;
             umap_res = project.comparisons.(comparison_name).dimension_reduction.umap;
         
             dim1 = pca_samp.score(umap_res.used_samples_idx,pc_x);
             dim2 = pca_samp.score(umap_res.used_samples_idx,pc_y);
             labels = "cluster " + string(project.comparisons.(comparison_name).dimension_reduction.kmeans.idx)';
            
             figure
             scatter(dim1,dim2,20,...
                     categorical(labels),'filled')
             add_labels(dim1, dim2, labels)
            
             title("PC "+ num2str(pc_x) +  " & " +  num2str(pc_y) + " with the cluster assignment from the kmeans!", 'FontSize',18)
             xlabel("PC" + num2str(pc_x) + " var: " + pca_samp.explained(pc_x), 'FontSize',18)
             ylabel("PC" + num2str(pc_y) + " var: " + pca_samp.explained(pc_y), 'FontSize',18)
             hold off
        else
             umap_res = project.comparisons.(comparison_name).dimension_reduction.umap;
             labels = project.comparisons.(comparison_name).dimension_reduction.kmeans.idx;
         
            dim1 = umap_res.reduction(:,1);
            dim2 = umap_res.reduction(:,2);
            labels = "cluster " + string(project.comparisons.(comparison_name).dimension_reduction.kmeans.idx)';
            dim_used = umap_res.used_pcs;
         
        
            figure
            scatter(dim1,dim2,20,...
                     categorical(labels),'filled')
            add_labels(dim1, dim2, labels)
            
            title("UMAP with the model lables", 'FontSize',18)
            xlabel("UMAP1 " + num2str(dim_used) + "PCs -> 2D", 'FontSize',14)
            ylabel("UMAP2 " + num2str(dim_used) + "PCs -> 2D", 'FontSize',14)
    
            hold off
        end


    end

    if ~isempty(rxn_to_visualize)

       
    end

end



function add_labels(dim1, dim2, labels)

    hold on
    
    unique_labels = unique(labels);
    data = [dim1 dim2];
    
    for i = unique_labels
        
        idx_label_samples = labels == i;
        pos_text = mean(data(idx_label_samples,:), 1);
        
        text(pos_text(1), pos_text(2), num2str(i), ...
            'VerticalAlignment','bottom', ...
            'HorizontalAlignment','right', ...
            'FontSize',18);
    end

end

