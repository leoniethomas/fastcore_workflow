function [kdl_matrix,p_value_kdl,sampling_sets,fdr] = perform_kdl_divergence_analysis(model,sampling_matrix,options)
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
        %   - sampling_matrix:  a sampling matrix to use as test dataset 
        %   - options:  
        %       - num_parallel_workers: workers you want to assign to this
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
        %       - number_of_ind_samplings: number of independent samplings,
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
            sampling_matrix =[]
            options.num_parallel_workers =10
            options.nPointsReturned =1000
            options.nStepsPerPoint =2000
            options.toRound =1
            options.number_of_ind_samplings =10
        end
        maxWorkers = parcluster('local').NumWorkers;
        disp("How many workers do I have ?")
        disp(maxWorkers);
        if maxWorkers-2 <= options.num_parallel_workers
            options.num_parallel_workers = maxWorkers -2;
        end
        
        sampling_sets = run_chrr_sampling(model,options,options.number_of_ind_samplings, options.num_parallel_workers);
        
        % distance between sets 
        pairs = nchoosek(1:options.number_of_ind_samplings, 2); 
        pairCell = num2cell(pairs, 2); % use as an input for arrayfun to run over all possible pairs between the n sets choosen
        pairwise_kdl = cell(length(pairs),1);

        delete(gcp('nocreate'))
        parpool(options.num_parallel_workers);
        parfor x=1:numel(pairCell)            
            pairwise_kdl{x} = get_kld_value_pairs(sampling_sets{pairCell{x}(1),1},sampling_sets{pairCell{x}(2),1});                          
        end
        pairwise_kdl = cell2mat(pairwise_kdl)'; 
        save('temp.mat'); % storing working env for debugging
        %load('temp.mat')
        
        if ~isempty(sampling_matrix)
            % distance between sets and sampling 

            pairCell_sampling_matrix = num2cell(options.number_of_ind_samplings, 2); % use as an input for arrayfun to run over all possible pairs between the n sets choosen
            pairwise_kdl_sampling_matrix = cell(options.number_of_ind_samplings,1);

            delete(gcp('nocreate'))
            parpool(options.num_parallel_workers);
            parfor x=1:numel(pairCell_sampling_matrix)            
                pairwise_kdl_sampling_matrix{x} = get_kld_value_pairs(sampling_sets{pairCell_sampling_matrix{x}(1),1},sampling_matrix);                          
            end
            pairwise_kdl_sampling_matrix = cell2mat(pairwise_kdl_sampling_matrix)'; 
            % pairwise_kdl_sampling_matrix = cellfun(@(x) get_kld_value_pairs(sampling_sets{x(1),1}, ...
            %                                                                 sampling_matrix),...
            %                                        pairCell_sampling_matrix, 'UniformOutput', false);
            % pairwise_kdl_sampling_matrix = cell2mat(pairwise_kdl_sampling_matrix)'; 
            train_data = pairwise_kdl;
            test_data = pairwise_kdl_sampling_matrix;
        else
            % in case sampling was not performed the disance measures of
            % the sets itself are split into train and test to get a
            % measure of how much we can trust a given rxn sampling
            % distribution

            N = size(pairwise_kdl, 2);        % total number of elements
            numSamples = round(0.1 * N);      % number of 1s (20%)
            binaryVec = zeros(1, N);
            randIdx = randperm(N, numSamples); % Randomly choose positions to set to 1
            binaryVec(randIdx) = 1;
            if sum(binaryVec) == 0
                binaryVec(1) = 1;
            end
            train_data = pairwise_kdl(:, find(~binaryVec));
            test_data = pairwise_kdl(:, find(binaryVec));
        end
        
        p_value_kdl = cell2mat(cellfun(@(rxn_idx) ranksum(train_data(rxn_idx(1),:),test_data(rxn_idx(1),:)),num2cell(1:size(pairwise_kdl, 1))',"UniformOutput",false));
        %p_adj_kdl = mafdr(p_value_kdl,'BHFDR', true);
        
        
        figure
        hist(p_value_kdl,100)
        sum( p_value_kdl < 0.05)
        fdr = sum( p_value_kdl < 0.05)/size(pairwise_kdl,1)
        if fdr > 0.05
            error("The False Discovery rate is lower than 5 % for the testing of sampling results that are obtained on the same model! Cause to worry, check again your samples, maybe increase sample number of samples sets?")
        end

        kdl_matrix = pairwise_kdl;



end

function [kld_vector] = get_kld_value_pairs(X,Y)

    if size(X,1) ~= size(Y,1)
        error("The two sampling results do not have the same dimensions, likely the rxns are then also not in the same order!!")
    end

    N = size(X,1);
    kld_vector = zeros(1, N);
    
    for rxn_idx = 1:N
        
        kld_vector(rxn_idx) = get_kld_value(X(rxn_idx,:), Y(rxn_idx,:));
    end
    % figure
    % histogram(kld_vector)

end

function [samples] = run_chrr_sampling(model,options,number_of_ind_samplings,num_parallel_workers);
    arguments
    model
    options
    number_of_ind_samplings
    num_parallel_workers =3
    end
    
    samples = cell(number_of_ind_samplings,1); 
    delete(gcp('nocreate'))
    parpool(num_parallel_workers); % 3 workers should be available on every laptop theoretically
    parfor x = 1:number_of_ind_samplings 
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



function kdl_value_rxn = get_kld_value(xdata,ydata)
    xmin = min(xdata);
    xmax = max(xdata);
    ymin = min(ydata);
    ymax = max(ydata);

    epsilon = 1e-12;  % tiny margin
    delta = 1e-6 * (xmax - xmin);   % small buffer
    Support = [xmin - delta - epsilon, xmax + delta + epsilon];
    
    [P.y, P.x] = ksdensity(xdata, ...
        'NumPoints', 5000, ...
        'Support', Support,...
        'BoundaryCorrection','reflection'); 

    delta = 1e-6 * (ymax - ymin);   % small buffer
    Support = [ymin - delta - epsilon, ymax + delta + epsilon];
    [Q.y, Q.x] = ksdensity(ydata, ...
        'NumPoints', 5000, ...
        'Support', Support,...
        'BoundaryCorrection','reflection'); 
   
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
    
    kdl_value_rxn = My_KLD(P, Q);

end

 function KLD_value = My_KLD(P, Q)
    % P and Q are structs with fields:
    %   P.x, P.y
    %   Q.x, Q.y
    
    p = P.y(:);
    q = Q.y(:);
    
    
    if all(q == 0)
        eps_val = min(p(p ~= 0)) * 0.1;
    else
        eps_val = min( ...
            min(p(p ~= 0)), ...
            min(q(q ~= 0)) ) * 0.1;
    end
    
    
    p(p < eps_val) = eps_val;
    q(q < eps_val) = eps_val;
    
    
    f = p .* log(p ./ q);
    
   
    dx = P.x(2) - P.x(1);
   
    KLD_value = dx * 0.5 * ...
        (f(1) + 2*sum(f(2:end-1)) + f(end));
end

