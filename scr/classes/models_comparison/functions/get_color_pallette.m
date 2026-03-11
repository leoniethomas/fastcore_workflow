function map = get_color_pallette()
        nColors = 256;        % number of color steps
        negColor = [0 0 1];   % blue
        midColor = [1 1 1];   % white
        posColor = [1 0 0];   % red
        
        % Interpolate negative to mid (blue → white)
        cmap_neg = [linspace(negColor(1), midColor(1), nColors/2)', ...
                    linspace(negColor(2), midColor(2), nColors/2)', ...
                    linspace(negColor(3), midColor(3), nColors/2)'];
        
        % Interpolate mid to positive (white → red)
        cmap_pos = [linspace(midColor(1), posColor(1), nColors/2)', ...
                    linspace(midColor(2), posColor(2), nColors/2)', ...
                    linspace(midColor(3), posColor(3), nColors/2)'];
        
        % Combine to single colormap
        cmap = [cmap_neg; cmap_pos];
        map = colormap(cmap)

end