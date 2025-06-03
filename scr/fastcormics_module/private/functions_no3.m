classdef functions_no3
   methods
      function rev_find = reverse_find_from_indices(obj,length,indices)
        rev_find = zeros(1,length)';
        rev_find(indices) = 1; 
      end

      function [final_counts] = get_pathway_counts(obj,AA,model_orig)
        subsystem_names =  unique([model_orig.subSystems{:}]);
        [pathways, ~, ub] = unique([model_orig.subSystems{find(AA),1}]);
        path_counts = histc(ub, 1:length(pathways));
        T = table(pathways', path_counts);
        [~, ia, ib] = intersect(subsystem_names, T.Var1);
        final_counts = zeros(1,length(subsystem_names))';
        final_counts(ia) = T.path_counts(ib);
        %final_counts = final_counts{:};
      end
      
      function [expanded_values] = expand_values(obj,idx,values,max_val)
          min_zero = zeros(max_val,1);
          min_zero(idx,1) = values;
          expanded_values = min_zero;
      end

      function [fig] = plot_clustergram(obj,data,rownames, colnames,title,position,altcolor)
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
      
   end
end