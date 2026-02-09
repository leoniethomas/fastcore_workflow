function visualize_sampling_landscape(project,comparison_name,pcs, rxn_to_visualize,sampling_feature,num_clusters)
    arguments
        project 
        comparison_name
        pcs (1,2) =[1,2]
        rxn_to_visualize (:,1) string =["biomass_reaction"]
        sampling_feature ="flux"
        num_clusters =0
    end
    if num_clusters == 0
        num_clusters = length(unique(project.comparisons.(comparison_name).sample_model_labels));
    end
    reference_model = project.comparisons.(comparison_name).reference_model;
    
    num_clusters_kmeans = num_clusters;
    sample_model_labels = project.comparisons.(comparison_name).sample_model_labels;
    if sampling_feature == "flux"
        ordered_samples = project.comparisons.(comparison_name).ordered_samples;
        dim_names = project.models.(reference_model).model.rxns;
    elseif sampling_feature == "fluxsum"
        ordered_samples = project.comparisons.(comparison_name).ordered_samples_fluxsum;
        dim_names = project.models.(reference_model).model.mets;
    else
        error("You choose a  non valid argument for the sampling feature! Choose flux or fluxsum, by default the flux per reaction will be portrayed!")
    end
    pc_x = pcs(1);
    pc_y = pcs(2);
    rxn_color = rxn_to_visualize;

    
    % perform a pca 
    disp("compute pcs!")
    pca_samp  = struct();
    % filter out dimensions which are always 0
    ordered_samples_filter = ordered_samples(find(any(ordered_samples ~= 0,2)),:);
    ordered_dim_names = dim_names(find(any(ordered_samples ~= 0,2)),:);
    [pca_samp.coeff,pca_samp.score,pca_samp.latent,pca_samp.tsquared,pca_samp.explained] = pca(ordered_samples_filter');
    disp("compute kmeans!")
    km = kmeans(ordered_samples_filter',num_clusters_kmeans);
    [s,~] = silhouette(ordered_samples_filter',km,'Euclidean');
    mean_sil = mean(s);
    %
    homogen = [];
    for k = unique(km)'
        labels_in_cluster = sample_model_labels(find(km == k));
        [s,~,j]=unique(labels_in_cluster);
        f = s{mode(j)};
        m = sum(labels_in_cluster == f);
        homogen = [ homogen, m/length(labels_in_cluster)];
    end
    homogen = mean(homogen);
    
    
    figure
    scatter(pca_samp.score(:,pc_x),pca_samp.score(:,pc_y),20,km,'filled')
    
    for i = unique(km)'
        disp(i)
        idx_label_samples = find(km == i);
        %disp(mean(idx_label_samples))
        pos_text = mean(pca_samp.score(idx_label_samples,:));
        text(pos_text(pc_x),pos_text(pc_y), num2str(i), 'VerticalAlignment',...
             'bottom', 'HorizontalAlignment', 'right', 'FontSize',18);
        hold on 
    end
    
    title("PC "+ num2str(pc_x) +  " & " +  num2str(pc_y) + " with the cluster ids from kmeans", 'FontSize',18)
    xlabel("PC" + num2str(pc_x) + " var: " + pca_samp.explained(pc_x), 'FontSize',18)
    ylabel("PC" + num2str(pc_y) + " var: " + pca_samp.explained(pc_y), 'FontSize',18)
    hold off
    
    figure
    scatter(pca_samp.score(:,pc_x),pca_samp.score(:,pc_y),20,categorical(sample_model_labels),'filled')
    
    unique_labels = categories(categorical(sample_model_labels));
    for i = 1:numel(unique_labels)
    
        l = string(unique_labels(i));
        idx_label_samples = find(sample_model_labels == l);
        disp(mean(idx_label_samples))
        pos_text = mean(pca_samp.score(idx_label_samples,:));
        text(pos_text(pc_x),pos_text(pc_y), unique_labels(i), ...
             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right',...
             'FontSize',18);
        hold on 
    end
    title("PC "+ num2str(pc_x) +  " & " +  num2str(pc_y) + " with the condition labels", 'FontSize',18)
    xlabel("PC" + num2str(pc_x) + " var: " + pca_samp.explained(pc_x), 'FontSize',18)
    ylabel("PC" + num2str(pc_y) + " var: " + pca_samp.explained(pc_y), 'FontSize',18)
    hold off
    
    
    for rxn_name = rxn_color'
        
        
        rxn_id = find(matches(ordered_dim_names,rxn_name));
        if isempty(rxn_id)
            display("This metabolite/reaction is not used in any of the samples:" + rxn_name)
            continue
        end
        figure
        scatter(pca_samp.score(:,pc_x),pca_samp.score(:,pc_y),30,...
                ordered_samples_filter(rxn_id,:),'filled',...
                'MarkerFaceAlpha',0.5)
        
        unique_labels = categories(categorical(sample_model_labels));
        for i = 1:numel(unique_labels)
        
            l = string(unique_labels(i));
            idx_label_samples = find(sample_model_labels == l);
            disp(mean(idx_label_samples));
            pos_text = mean(pca_samp.score(idx_label_samples,:));
            text(pos_text(pc_x),pos_text(pc_y), unique_labels(i), ...
                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right',...
                 'FontSize',18);
            hold on 
        end
        title("PC "+ num2str(pc_x) +  " & " +  num2str(pc_y) + " with the condition labels",'FontSize',18)
        xlabel("PC" + num2str(pc_x) + " var: " + pca_samp.explained(pc_x),'FontSize',18)
        ylabel("PC" + num2str(pc_y) + " var: " + pca_samp.explained(pc_y),'FontSize',18)
        set(gca,'FontSize',18)
        
        cb = colorbar;
        cb.Label.String = regexprep(rxn_name,"_", " ");
        cb.FontSize = 18;
        
        hold off
    end

end

