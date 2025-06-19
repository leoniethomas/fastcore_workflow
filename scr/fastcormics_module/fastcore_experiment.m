classdef fastcore_experiment
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sample_labels
        rxn_names
        samples
        run_names
        fastcore_runs
        transformed_samples
        fluxsum
        met_names
    end
    
    methods
        %%
        function obj = fastcore_experiment(sampling_files,run_fluxsum)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            arguments
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
        %%
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
        function obj = join_fluxsum_output(obj)
            all_mets = cellfun(@(x) obj.fastcore_runs.(x).model.mets, string(fieldnames(obj.fastcore_runs)),'UniformOutput',false);
            all_mets = unique(vertcat(all_mets{:}));
            
            x = 'samplingResults_MDA_MB231_Cont_NO_model_20250602_090252';
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
        
    end
end

function mustBeMember_inModel(models_to_compare,obj)
            assert(all(ismember(models_to_compare,obj.run_names)), "The specified models to compare are not part of the model! Check the obj.run_names slot to get the correct labels!")
end



function mustBeFileType(file_path,needed_file_format_ending)
            assert(~isempty(regexp(file_path,needed_file_format_ending + "$")),...
                   "Input must be a " + needed_file_format_ending + " file!!")
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








