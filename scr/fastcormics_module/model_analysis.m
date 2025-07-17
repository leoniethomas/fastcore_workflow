classdef model_analysis
 
    properties
        model_names
        model_size
        reaction_presence
    end
    
    methods
        function obj = model_analysis(exp)
            
            
            condition_models = exp.condition_models;
            obj.model_names = fieldnames(condition_models);
            condition_models = structfun(@(x) removeUnusedGenesFastbox(x,1), ... % remove unused genes to get model size in terms of genes active 
                             condition_models,'UniformOutput',false);
                         
            % simple quantities per model - #rxns #genes #metabolites in the models
            obj.model_size = array2table(struct2array(structfun(@(x) {numel(x.rxns);numel(x.mets);numel(x.genes)}, ...
                                                                    condition_models,'UniformOutput',false))',...
                                             'VariableNames',{'count_reactions','count_metabolites','count_genes'},...
                                             'RowNames',obj.model_names);
            % get reaction precense
            AA_keep = struct2array(structfun(@(x) get_feature_presence(length(exp.original_model.rxns),x.AA), ...
                                   condition_models,'UniformOutput',false));
            obj.reaction_presence = AA_keep;
            
        end
        
        function fig = get_jaccard_similarity(obj)
           J = squareform(pdist(obj.reaction_presence','jaccard'));
           disp("Jaccard similarity:")
           1-J
           fig = plot_clustergram(1-J,...
                     obj.model_names,...
                     obj.model_names,...
                     {'Model similarity based on Jaccard distance of rxns existence in the model!'},...
                     [100 100 800 600]);   
        end
        
        function [fig,intersection,outersection] = get_intersections_outersections(exp,model_names,slot)
            arguments
               exp
               model_names 
               slot (1,1) string ="rxns"
            end
            
            exp.condition_models 
            
            
            
        end
    end
end

function [fig] = plot_clustergram(data,rownames, colnames,title,position,altcolor)
    arguments
        data
        rownames
        colnames
        title
        position
        altcolor =[255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; 
    end
  cgo_J = clustergram(data,...
              'RowLabels', rownames,...
              'ColumnLabels', colnames,...
              'ColumnLabelsRotate',45, ...
              'Cluster', 'all', ...
              'symmetric','False',...
              'Colormap', altcolor);  
   addTitle(cgo_J,title)
   cgf = plot(cgo_J); % This should be a figure handle
   colorbar(cgf,'eastoutside');
   fig = gcf;
   fig.Position = position;

end

 function rev_find = get_feature_presence(length,indices)
            rev_find = zeros(1,length)';
            rev_find(indices) = 1; 
 end



