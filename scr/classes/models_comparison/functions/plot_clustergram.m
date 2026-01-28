function [fig] = plot_clustergram(data,rownames, colnames,title,colorbarLabel,altcolor)
           
          cgo_J = clustergram(data,...
                      'RowLabels', rownames',...
                      'ColumnLabels', colnames',...
                      'ColumnLabelsRotate',45, ...
                      'Cluster', 'all', ...
                      'symmetric','False',...
                      'Standardize', 'none',...
                      'Colormap', altcolor);  
           addTitle(cgo_J,title)
           cgf = plot(cgo_J); % This should be a figure handle
           
           colorbar(cgf,'eastoutside');

           % Add colorbar and label
           cb = colorbar(cgf, 'eastoutside');
           cb.Label.String = colorbarLabel;
           % Set colorbar font size
           cb.Label.FontSize = 16;
           cb.FontSize = 16;  % tick labels        

           fig = gcf;
           
           ax = findall(fig, 'Type', 'Axes');
           set(ax, 'FontSize', 16);
           txt = findall(fig, 'Type', 'Text');
           set(txt, 'FontSize', 16)

        end
    
       