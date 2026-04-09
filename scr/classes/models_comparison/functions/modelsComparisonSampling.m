function project = modelsComparisonSampling(project,comparison_name)


    list_model_names = strsplit(comparison_name, "__"); 
    list_model_names = strsplit(list_model_names(1),"_vs_");
    models_list = rmfield(project.models, setdiff(fieldnames(project.models), list_model_names));
    % give the comparison the name of all models + a identifier choosen
    reference_model = project.comparisons.(comparison_name).reference_model;
    
    % run structural model comparison
    replacement_value = "analysis.sampling.samples"; % get the fba solution values
    project.comparisons.(comparison_name).ordered_samples = getOrderedFeatureMatrix(project,list_model_names,"rxns",reference_model,replacement_value);
    replacement_value = "analysis.FBA.v"; % get the fba solution values
    project.comparisons.(comparison_name).ordered_fba = getOrderedFeatureMatrix(project,list_model_names,"rxns",reference_model,replacement_value);
    
    sample_count_models = structfun(@(x) size(x.analysis.sampling.samples,2), models_list);
    project.comparisons.(comparison_name).sample_model_labels = repelem(list_model_names, sample_count_models);
    project.comparisons.(comparison_name).plots.sampling = visualize_sampling_landscape(project,comparison_name,'visible_plot',"off")


end


    




%% old wilcoxon test
% function stats = diff_flux_testing(samples_model1, samples_model2, reference_model)
% 
% 
%             stats = zeros(length(reference_model.model.rxns), 1);  % Preallocate for speed
%             for counter=1:length(reference_model.model.rxns)
%                 rxn_sample_value_model1 = samples_model1(counter,:);
%                 rxn_sample_value_model2 = samples_model2(counter,:);    
%                 %stats(counter) = ranksum(rxn_sample_value_model1,rxn_sample_value_model2 );
%                 %[h,stats(counter),ksstat] = kstest2(rxn_sample_value_model1,rxn_sample_value_model2);
%                 stats(counter) = ttest2(rxn_sample_value_model1,rxn_sample_value_model2 );
%                 A = rxn_sample_value_model2;
%                 B = rxn_sample_value_model1;
%                 direction_change_m = zeros(10600, 1);
%                 if (A >= 0 & B >= 0) | (A <= 0 & B <= 0)
%                     % Mark a change in direction (difference is zero)
%                     direction_change_m(counter) = 0;
%                 elseif A <= 0 & B >= 0
%                     % Mark a change in direction
%                     direction_changem(counter) = 1;
%                     % snr =(mean(B) - mean(A)) / (std(A) + std(B)); % signal to noise ratio
%                     % snr =abs(snr);
%                 elseif A >= 0 & B <= 0
%                     % Mark a change in direction
%                     direction_changem(counter) = -1;
%                     % snr =(mean(B) - mean(A)) / (std(A) + std(B)); % signal to noise ratio
%                     % snr = abs(snr);
%                 else
%                     disp('hello')
%                     % hist(A,100)
%                     % hist(B,100)
%                 end
%             end
% 
%             clean_idx = find(~isnan(stats));
% 
%             stats = stats(clean_idx);
%             [p_adj] = mafdr(stats,'BHFDR', true);
%             stats = [stats,p_adj];
% 
%             samples_model1 = samples_model1(clean_idx,:);
%             samples_model2 = samples_model2(clean_idx,:);
% 
%             % compute signal to noise ration 
%             % + adding 1000 to all samples flux values, so that we are only dealing on
%             % a positive scale 
%             % + we are doing the snr ratio -> deviding through the std of both
%             % distributions -> normally only the noise std 
%             % + but here we do not really have a baseline and a noise distribution 
%             mean_model1 = mean(samples_model1 + 1000, 2)-1000;  % mean of each row (observation)
%             median_model1 = median(samples_model1 + 1000, 2)-1000;  % median of each row (observation)
%             std_model1 = std(samples_model1 + 1000, 0, 2);  % std of each row (observation)
% 
%             mean_model2 = mean(samples_model2 + 1000, 2)-1000;  % mean of each row (observation)
%             median_model2 = median(samples_model2 + 1000, 2)-1000;  % mean of each row (observation)
%             std_model2 = std(samples_model2 + 1000, 0, 2);  % std of each row (observation)
% 
%             snr = (mean_model2 - mean_model1) ./ (std_model1 + std_model2);
% 
%             log2FC=log2(abs(mean_model2./mean_model1));
%             stats =[mean_model1, mean_model2,median_model1, median_model2, std_model1, std_model2, ...
%                           mean_model1 - mean_model2, log2FC, snr, ...
%                           stats];
%             stats = array2table(stats,"RowNames",reference_model.model.rxns(clean_idx),...
%                 "VariableNames",["mean model1", "mean model2", "median model1", "median model2", ...
%                  "std_model1", "std_model2", "mean_model1 - mean_model2", "log2FC", "snr","pval", "fdr"]);
%             stats = stats(find(~isnan(stats.log2FC)),:);
% 
% 
% 
%             % if writetofile
%             %     stats_tab=array2table(stats,'RowNames',obj.rxn_names,'VariableNames',{'mean_model1','mean_model2','median_model1','median_model2','std_model1','std_model2','diff','log2FC','SNR','pValue','p_adj'});
%             %     writetable(stats_tab, filepath + filename,'WriteRowNames',true)
%             % end
%             % if figflag
%             %     figure
%             %     histogram(-log2(stats(:,11)))
%             %     %title('P values (-log10)')
%             %     title( 'Histogram P adj values (-log2)')
%             % 
%             %     figure
%             %     %hist(log10(stats(:,4)))
%             %     histogram(stats(:,8))
%             %     %title('log10 foldchange (mean(B)/mean(A))')
%             %     title('Histogram log2 foldchange mean\_model2/mean\_model1))')
%             % 
%             %     figure
%             %     % plot(log10(stats(:,4)),stats(:,6),'*')
%             %     % title('vulcano: log10 foldchange vs -log10(P)')
%             %     plot(stats(:,8),-log2(stats(:,11)),'.')
%             %     hold on 
%             %     yline(-log2(0.001))
%             %     xline(1.5)
%             %     xline(-1.5)
%             %     text(-10,100,["adj\_pval < 0.001", "abs(logFC) > 1.5"])
%             %     title('vulcano: log2 foldchange vs -log2(P\_adj)')
%             %     xlabel("log2FC flux")
%             %     ylabel("-log2(p\_adj)")
%             % end
% 
%         end