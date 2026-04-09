function newFig = show_figure(figHandle)
% DUPLICATEFIGURE Creates a copy of a figure safely.
%   newFig = duplicateFigure(figHandle) duplicates the figure given by
%   figHandle. Works for figures with UIAxes, tables, clustergrams, etc.
%
% Example:
%   newFig = duplicateFigure(plots.funct.import);
    try
        
        % tf = any(arrayfun(@(x) isa(x, 'matlab.ui.control.UIAxes'), figHandle.Children));
        % if tf 
        %     error("This is a ui figure it needs to be saved and reopened!")
        % end
        figure(copyobj(figHandle, groot));
        return;  % success, exit function
   
    catch
        try
            % if tf 
            %     error("This is a ui figure it needs to be saved and reopened!")
            % end
            children = get(figHandle, 'Children');  % get figure children
            newFig = figure;                        % create new figure
            copyobj(children, newFig);              % attempt copy
            set(newFig, 'Visible', 'on');           % ensure it shows
            return;  % success, exit function
        catch 
            % try
            %     % if tf 
            %     %     error("This is a ui figure it needs to be saved and reopened!")
            %     % end
                % Check input
                if ~ishandle(figHandle) || ~strcmp(get(figHandle,'Type'),'figure')
                    error('Input must be a valid figure handle.');
                end
                
            
                % Generate a temporary filename in the system temp folder
                tempFile = [tempname, '.fig'];
            
                try
                    % Save the figure to the temporary file
                    savefig(figHandle, tempFile);
            
                    % Open the figure as a new figure
                    newFig = openfig(tempFile, 'reuse'); % opens a new figure window
            
                    % Ensure it is visible
                    set(newFig, 'Visible', 'on');
            
                catch ME
                    % If any error occurs, delete the temp file and rethrow
                    if exist(tempFile, 'file')
                        delete(tempFile);
                    end
                    rethrow(ME);
                end
            
                % Delete the temporary file
                if exist(tempFile, 'file')
                    delete(tempFile);
                end
            % catch
            %     drawnow
            %     frame = getframe(figHandle);
            %     figure('Position', [10 10 1000 1000]) % [left bottom width height]
            %     imshow(frame.cdata)
            % end
        end
    end

end