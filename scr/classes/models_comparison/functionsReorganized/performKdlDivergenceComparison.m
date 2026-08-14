function [kdl_matrix,p_value_kdl,sampling_sets,fdr] = performKdlDivergenceComparison(model,samplingMatrix,options)
        % This function performs an evaluation of the variability and
        % convergence of the estimated sampling distribution. This
        % estimation is meant to give a measure of how trustworthy the
        % sampling results are and therefore should lower the False
        % Discovery rate. The method was adpoted from Bruno G. Galuzzi et al 2024 
        %-> link to paper: https://www.sciencedirect.com/science/article/pii/S1532046424000157?via%3Dihub
        % and the code was written in accordance to the analysis described
        % in the sections 2.5 and 3.4 of the publication. 
        % The basic question this function asks is which of the
        % distributions that we observed due to sampling can be trusted ? 
        % are there rxns for which the sampling has not converged. 
        % This is relevant since the rxn distribution needs to be converged
        % in order for a downstream differential analysis between
        % conditions to be trustworthy. The provided pvalues can serve as a
        % measure to filter for "trustworthy" rxns. 
        % The aim is to push the False discovery rate below 5 percent with
        % this method.
        % The basis of the method is the Kullback-Leibler divergence which
        % measures the difference between two probability distributions. 
        % In the function below multiple sampling rounds are performed and
        % then the pairwise KLD are computed. The distances are splitted
        % into train and test and then a ranksum test is performed to see
        % wether for all of the rxns we get a pvalue of above 0.05. If this
        % is not the case then the test indicates that the distributions of
        % samplings done on the same model are significanly different,
        % which should not be the case and is an indice for the non
        % convergence of this specific rxns sampling distribution. 
        % Input: 
        %   - model:            cobra model object 
        %   - samplingMatrix:  a sampling matrix to use as test dataset 
        %   - options:  
        %       - numParallelWorkers: workers you want to assign to this
        %                               part in the parfor loop (default is
        %                               10 if the computer has at least 12
        %                               workers), otherwise set to the
        %                               amount of workers available -2 
        %       - nPointsReturned:      how many samples to draw (default
        %                               is 1000)
        %       - nStepsPerPoint:       thinning parameter, how many steps
        %                               do I have to make in my markov chain
        %                               hit and run to draw one sample 
        %                               (default: draw every 2000th step as
        %                               a sample). The higher this value
        %                               the lower the autocorrelation
        %                               between the samples, theoretically
        %       - toRound:              set to one for rounding the
        %                               solution space before performing the
        %                               hit-and-run, leads in theory to the
        %                               hit-and-run not to be trapped in the corners
        %                               and elongated forms of the flux polytope
        %       - numberOfIndSamplings: number of independent samplings,
        %       the more you have the better the estimation of your
        %       pairwise KLdivergence, the more trust you can put into the
        %       pvalue at the end of this function (default: 10, should be
        %       sufficient, leads to 45 pairsise KLD metrics per rxns), in
        %       the original publication 20 sets were generated)
        %
        % Output: 
        %   - kdl_matrix:   The pairwise kld metric per rxns (#rxns x #
        %                   pairwise set comparisons)
        %   - p_value:      measeure of how likely it is that the sampling
        %                   can be trusted, ranksum pvalue which gives you
        %                   the likelyhood of the samplings not coming from
        %                   the same solution space (# rxns x 1)
        %   - sampling_sets:The sampling sets drawn from the solution space
        %                   cell array, every cell array entails a matrix with
        %                   dims (#rxns x # samples)
        %                   
        
        arguments
            model
            samplingMatrix =[]
            options.numParallelWorkers =10
            options.nPointsReturned =1000
            options.nStepsPerPoint =2000
            options.toRound =1
            options.numberOfIndSamplings =10
        end
        maxWorkers = parcluster('local').NumWorkers;
        disp("How many workers do I have ?")
        disp(maxWorkers);
        if maxWorkers-2 <= options.numParallelWorkers
            options.numParallelWorkers = maxWorkers -2;
        end
        
        sampling_sets = runChrrSampling(model,options,options.numberOfIndSamplings, options.numParallelWorkers);
        
        % distance between sets 
        pairs = nchoosek(1:options.numberOfIndSamplings, 2); 
        pairCell = num2cell(pairs, 2); % use as an input for arrayfun to run over all possible pairs between the n sets choosen
        pairwiseKdl = cell(length(pairs),1);

        delete(gcp('nocreate'))
        parpool(options.numParallelWorkers);
        parfor x=1:numel(pairCell)            
            pairwiseKdl{x} = getKldValuePairs(sampling_sets{pairCell{x}(1),1},sampling_sets{pairCell{x}(2),1});                          
        end
        pairwiseKdl = cell2mat(pairwiseKdl)'; 
        save('temp.mat'); % storing working env for debugging
        %load('temp.mat')
        
        if ~isempty(samplingMatrix)
            % distance between sets and sampling 
            pairCell_samplingMatrix = num2cell((1:options.numberOfIndSamplings)');
            pairwiseKdl_samplingMatrix = cell(options.numberOfIndSamplings,1);

            delete(gcp('nocreate'))
            parpool(options.numParallelWorkers);
            parfor x=1:numel(pairCell_samplingMatrix)            
                pairwiseKdl_samplingMatrix{x} = getKldValuePairs(sampling_sets{pairCell_samplingMatrix{x}(1),1},samplingMatrix);                          
            end
            pairwiseKdl_samplingMatrix = cell2mat(pairwiseKdl_samplingMatrix)'; 
            % pairwiseKdl_samplingMatrix = cellfun(@(x) getKldValuePairs(sampling_sets{x(1),1}, ...
            %                                                                 samplingMatrix),...
            %                                        pairCell_samplingMatrix, 'UniformOutput', false);
            % pairwiseKdl_samplingMatrix = cell2mat(pairwiseKdl_samplingMatrix)'; 
            train_data = pairwiseKdl;
            test_data = pairwiseKdl_samplingMatrix;
        else
            % in case sampling was not performed the disance measures of
            % the sets itself are split into train and test to get a
            % measure of how much we can trust a given rxn sampling
            % distribution

            N = size(pairwiseKdl, 2);        % total number of elements
            numSamples = round(0.3 * N);      % number of 1s (20%)
            binaryVec = zeros(1, N);
            randIdx = randperm(N, numSamples); % Randomly choose positions to set to 1
            binaryVec(randIdx) = 1;
            if sum(binaryVec) == 0
                binaryVec(1) = 1;
            end
            train_data = pairwiseKdl(:, find(~binaryVec));
            test_data = pairwiseKdl(:, find(binaryVec));
        end
        
        p_value_kdl = cell2mat(cellfun(@(rxn_idx) ranksum(train_data(rxn_idx(1),:),test_data(rxn_idx(1),:)),num2cell(1:size(pairwiseKdl, 1))',"UniformOutput",false));
        %p_adj_kdl = mafdr(p_value_kdl,'BHFDR', true);
        
        
        figure
        histogram(p_value_kdl, 100)
        sum( p_value_kdl < 0.05)
        fdr = sum( p_value_kdl < 0.05)/size(pairwiseKdl,1)
        % if fdr > 0.05
        %     error("FDR is %.1f%% — above the 5%% threshold. Sampling may not have converged.", fdr*100)
        % end

        kdl_matrix = pairwiseKdl;



end

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

end

function [samples] = runChrrSampling(model,options,numberOfIndSamplings,numParallelWorkers);
    arguments
    model
    options
    numberOfIndSamplings
    numParallelWorkers =3
    end
    
    samples = cell(numberOfIndSamplings,1); 
    delete(gcp('nocreate'))
    parpool(numParallelWorkers); % 3 workers should be available on every laptop theoretically
    parfor x = 1:numberOfIndSamplings 
        pause(30*x); % permits the different loops to go into the cobra directory at the same time and init, 
        % doing it at the exact parallel causes issues & the script to crash
        initCobraToolbox(false) 
        changeCobraSolver("gurobi")
        
        %%%% IMPORTANT, never remove the random seed from the sampling!!!
        rng(randi([1, 1000])); % random seed to get different samples
        %%%% --------------------------
        sampleFile =  char(string(datetime("now", "Format", "yyyyMMdd_HHmm")));
        [~, samples{x}] = sampleCbModel(model, sampleFile,  'CHRR', options);
    end
end



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

