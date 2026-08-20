function [kld_vector] = getKldValuePairs(X,Y)

    if size(X,1) ~= size(Y,1)
        error("The two sampling results do not have the same dimensions, likely the rxns are then also not in the same order!!")
    end

    N = size(X,1);
    kld_vector = zeros(1, N);
    
    for rxn_idx = 1:N
        kld_vector(rxn_idx) = getKLDValue(X(rxn_idx,:), Y(rxn_idx,:));
    end
    % figure
    % histogram(kld_vector)


    function kdlValueRxn = getKLDValue(xdata,ydata)
        % Add at the top of getKLDValue:
        if var(xdata) == 0 && var(ydata) == 0
            kdlValueRxn = 0;  % identical constant distributions
            return
        end
        if var(xdata) == 0 || var(ydata) == 0
            kdlValueRxn = NaN;  % undefined — flag for downstream filtering
            return
        end
        xdata = double(xdata);
        ydata = double(ydata);
        allData = [xdata, ydata];
        xmin = min(allData);
        xmax = max(allData);
        epsilon = 1e-12;
        delta   = 1e-6 * (xmax - xmin);
        sharedSupport = [xmin - delta - epsilon, xmax + delta + epsilon];
        sharedGrid    = linspace(sharedSupport(1), sharedSupport(2), 5000);
    
        P.y = ksdensity(xdata, sharedGrid, 'Support', sharedSupport, ...
                        'NumPoints', 5000, ...
                        'BoundaryCorrection', 'reflection');
        Q.y = ksdensity(ydata, sharedGrid, 'Support', sharedSupport, ...
                        'NumPoints', 5000, ...
                        'BoundaryCorrection', 'reflection');
       
        % figure;
        % 
        % % ---- Subplot 1: xdata ----
        % subplot(2,1,1);  % 2 rows, 1 column, first subplot
        % histogram(xdata, 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8], 'NumBins', 100);
        % hold on;
        % plot(P.x, P.y, 'r-', 'LineWidth', 2);
        % xlabel('xdata values');
        % ylabel('Probability Density');
        % title('Histogram and KDE of xdata');
        % legend('Histogram', 'KDE');
        % grid on;
        % 
        % % ---- Subplot 2: ydata ----
        % subplot(2,1,2);  % 2 rows, 1 column, second subplot
        % histogram(ydata, 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8], 'NumBins', 100); 
        % hold on;
        % plot(Q.x, Q.y, 'r-', 'LineWidth', 2);
        % xlabel('ydata values');
        % ylabel('Probability Density');
        % title('Histogram and KDE of ydata');
        % legend('Histogram', 'KDE');
        % grid on;
    
        kdlValueRxn = KLDis(P.y, Q.y); 
    
    end


end


