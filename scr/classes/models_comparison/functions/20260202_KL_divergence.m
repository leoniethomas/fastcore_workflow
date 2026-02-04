%% Performing the FDR correction shown by Bruno G. Galuzzi et al 2024 
%-> link to paper: https://www.sciencedirect.com/science/article/pii/S1532046424000157?via%3Dihub


initCobraToolbox(false) 
changeCobraSolver("gurobi")

working_path = "/Users/leonie.thomas/Documents/fastcore_workflow_with_vanille";
cd (working_path)
addpath(genpath(working_path))

load(working_path + filesep + "context_specific_models" + filesep + "20260119_1042" + filesep + "project_23012026_1453_28012026_1508_obj_vanille_sampling.mat")
model = project.models.WT.model;

[kdl_matrix,p_value_kdl] = perform_kdl_divergence_analysis(model);





%% -------------------- DEFINE FUNCTIONS -------------------


function [kdl_matrix,p_value_kdl,sampling_sets] = perform_kdl_divergence_analysis(model,sampling_matrix,options)
        % This function 
        arguments
            model
            sampling_matrix =[]
            options.num_parallel_workers =10
            options.nPointsReturned =3000
            options.nStepsPerPoint =2000
            options.toRound =1
            options.number_of_ind_samplings =10
        end
        maxWorkers = parcluster('local').NumWorkers;
        disp("How many workers do I have ?")
        disp(maxWorkers);
        if maxWorkers < options.num_parallel_workers
            options.num_parallel_workers = maxWorkers -2;
        end
        
        sampling_sets = run_chrr_sampling(model,options,options.number_of_ind_samplings, options.num_parallel_workers);

        % distance between sets 
        pairs = nchoosek(1:options.number_of_ind_samplings, 2); 
        pairCell = num2cell(pairs, 2); % use as an input for arrayfun to run over all possible pairs between the n sets choosen
        pairwise_kdl = cellfun(@(x) get_kld_value_pairs(sampling_sets{x(1),1}, ...
                                                        sampling_sets{x(2),1}),...
                               pairCell, 'UniformOutput', false);
        pairwise_kdl = cell2mat(pairwise_kdl)'; 
        
        if ~isempty(sampling_matrix)
            % distance between sets and sampling 

            pairCell_sampling_matrix = num2cell(options.number_of_ind_samplings, 2); % use as an input for arrayfun to run over all possible pairs between the n sets choosen
            pairwise_kdl_sampling_matrix = cellfun(@(x) get_kld_value_pairs(sampling_sets{x(1),1}, ...
                                                                            sampling_matrix),...
                                                   pairCell_sampling_matrix, 'UniformOutput', false);
            pairwise_kdl_sampling_matrix = cell2mat(pairwise_kdl_sampling_matrix)'; 
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
            train_data = pairwise_kdl(:, find(~binaryVec));
            test_data = pairwise_kdl(:, find(binaryVec));
        end
        
        p_value_kdl = cell2mat(cellfun(@(rxn_idx) ranksum(train_data(rxn_idx(1),:),test_data(rxn_idx(1),:)),num2cell(1:size(pairwise_kdl, 1))',"UniformOutput",false));
        %p_adj_kdl = mafdr(p_value_kdl,'BHFDR', true);
        
        
        figure
        hist(p_value_kdl,100)
        sum( p_value_kdl < 0.05)
        FDR = sum( p_value_kdl < 0.05)/size(pairwise_kdl,1)
        if FDR > 0.05
            error("The False Discovery rate is lower than 5 % for the testing of sampling results that are obtained on the same model! Cause to worry, check again your samples, maybe increase sample number of samples sets?")
        end

        kdl_matrix = pairwise_kdl;



end

function [kld_vector] = get_kld_value_pairs(X,Y)

    if any(size(X) ~= size(Y))
        error("The two sampling results do not have the same dimensions, likely the rxns are then also not in the same order!!")
    end

    N = size(X,1);
    kld_vector = zeros(1, N);
    
    parfor rxn_idx = 1:N
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

        rng(randi([1, 1000])); % random seed to get different samples
        sampleFile =  char(string(datetime("now", "Format", "yyyyMMdd_HHmm")));
        [~, samples{x}] = sampleCbModel(model, sampleFile,  'CHRR', options);
    end
end



function kdl_value_rxn = get_kld_value(xdata,ydata)
    xmin = min(xdata);
    xmax = max(xdata);
    ymin = min(ydata);
    ymax = max(ydata);
    delta = 1e-6 * (xmax - xmin);   % small buffer
    [P.y, P.x] = ksdensity(xdata, ...
        'NumPoints', 5000, ...
        'Support', [xmin - delta, xmax + delta],...
        'BoundaryCorrection','reflection'); 
    [Q.y, Q.x] = ksdensity(ydata, ...
        'NumPoints', 5000, ...
        'Support', [ymin - delta, ymax + delta],...
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

