    function [mean_sil,homogen] = visualize_sampling(model_names, model_labels,samples, num_clusters_kmeans, ...
                                 pc_x, pc_y,disp_fig)
%     model_names = string(fieldnames(models))';
%     model_labels = regexprep(regexprep(string(fieldnames(condition_models)),"_"," "),...
%                              "MDA MB231 ","");
%     samples = samples_joined_ordered';
%     num_clusters_kmeans = numel(model_names);
%     pc_x = 2;
%     pc_y = 3;

    %%

    sample_labels = arrayfun(@(x) repmat(x,1,2000),model_names, ...
                             'UniformOutput',false);
    model_labels_all =  arrayfun(@(x) repmat(x,1,2000),model_labels, ...
                             'UniformOutput',false);                     
    sample_labels = [sample_labels{:}];
    model_labels_all = [model_labels_all{:}];
   
%     biomass_idx = 3625;
%     
%     sample_labels = sample_labels(find(samples(:,biomass_idx) > 0.25));
%     model_labels_all = model_labels_all(find(samples(:,biomass_idx) > 0.25));
%     samples = samples(find(samples(:,biomass_idx) > 0.25),:);
    %%
    [coeff,score,latent,tsquared,explained] = pca(samples);
    km = kmeans(samples,num_clusters_kmeans);
    [s,~] = silhouette(samples,km,'Euclidean');
    mean_sil = mean(s);
    %%
    homogen = [];
    for k = unique(km)'
        labels_in_cluster = sample_labels(find(km == k));
        [s,~,j]=unique(labels_in_cluster);
        f = s{mode(j)};
        m = sum(labels_in_cluster == f);
        homogen = [ homogen, m/length(labels_in_cluster)];
    end
    homogen = mean(homogen);

    %%
    if disp_fig
        figure
        hist(categorical(model_labels_all' + "__" + km))
        title("Count of samples per condition - kmeans cluster id combination")
        xlabel("condition__cluster_id")
        ylabel("count of samples")
        xtickangle(45)  
        %%
        figure
        %id_biomass_ordered = 3625;
        %scatter(score(:,pc_x),score(:,pc_y),5,samples(:,id_biomass_ordered))
        scatter(score(:,pc_x),score(:,pc_y),5,categorical(sample_labels))

        unique_labels = categories(categorical(sample_labels));
        for i = 1:numel(unique_labels)

            l = string(unique_labels(i));
            idx_label_samples = find(sample_labels == l);
            %disp(mean(idx_label_samples))
            pos_text = mean(score(idx_label_samples,:));
            text(pos_text(pc_x),pos_text(pc_y), model_labels(i), ...
                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
            hold on 
        end
        title("PC "+ num2str(pc_x) +  " & " +  num2str(pc_y) + " with the condition labels")
        xlabel("PC" + num2str(pc_x) + " var: " + explained(pc_x))
        ylabel("PC" + num2str(pc_y) + " var: " + explained(pc_y))
        hold off
    end

    %%
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

        title("PC "+ num2str(pc_x) +  " & " +  num2str(pc_y) + " with the cluster ids from kmeans")
        xlabel("PC" + num2str(pc_x) + " var: " + explained(pc_x))
        ylabel("PC" + num2str(pc_y) + " var: " + explained(pc_y))
        hold off
    end
end