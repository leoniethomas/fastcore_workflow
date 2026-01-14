classdef model_analysis
    % This object class is meant to store the analysis results that are
    % obtained when analysing and comparing context specific metabolic
    % models. 
 
    properties
        model_names % names of context specific models, for example named after SampleID or Treatment regime
        pathway_counts
        gene_essentiality
        enrichment
        fastcore_models
        fba
        fva
    end
    
    methods
        function obj = model_analysis(exp)
            % this function performs multiple tasks
            % - reading in context specific models, with names 
            % - getting the sice of the networks (genes,reaction, &
            % metabolite count)
            % - aligning the retainedRxns/retained_reaction output coming from
            % fastcormics function, and defining a new matrix which gives
            % the reacion presence of the models for all the reaction that
            % were in the original model
            arguments
               exp (1,1) struct
            end
            
            obj.fastcore_models = exp;
            obj.model_names = fieldnames(obj.fastcore_models.condition_models);
            
        end

        
        
        function disp(obj)
            % silently process all condition models
            condition_models = structfun(@(x) suppressOutput(x), ...
                                         obj.fastcore_models.condition_models, ...
                                         'UniformOutput', false);
            
            % convert to table
            data = struct2array(structfun(@(x) {numel(x.rxns); numel(x.mets); numel(x.genes)}, ...
                                          condition_models, 'UniformOutput', false))';
            
            T = array2table(data, ...
                            'VariableNames', {'count_reactions','count_metabolites','count_genes'}, ...
                            'RowNames', obj.model_names);
            disp(T)
        
            % --- nested function to suppress output ---
            function model = suppressOutput(x)
                [~, model] = evalc('removeUnusedGenesFastbox(x,1)');
            end
        end

        function react_presence = reaction_presence(obj)
            % get reaction precense
            react_presence = struct2array(structfun(@(x) get_feature_presence(length(obj.fastcore_models.original_model.rxns),x.retainedRxns), ...
                                   obj.fastcore_models.condition_models,'UniformOutput',false));
        end

        function enrichment = perform_gene_enrichment(obj,exp, threshold)
            %threshold = 0.5;
            data = obj.gene_essentiality.ratio;
            mask = data < threshold;

            num_cols = size(data, 2);
            row_indices = cell(1, num_cols);  % store indices for each column
            enrichment = struct()
            
            [found,idx] = ismember(regexprep(exp.original_model.genes, '\.\d+$', ''),exp.dico.ENTREZ);
            GeneList_all =  exp.dico.SYMBOL(idx(found),:);

            for col = 1:num_cols
                row_indices{col} = exp.original_model.genes(find(mask(:, col)));  % find row indices where condition is true
                [~,idx] = ismember(regexprep(row_indices{col}, '\.\d+$', ''),exp.dico.ENTREZ);
                row_indices{col} = exp.dico.SYMBOL(idx,:);
                enrichment.(obj.model_names{col}) = GeneEnrichments(string(row_indices{col}),GeneList_all);
            end
            
        end
        
        function [fig,J] = jaccard_similarity(obj,slot)
           % this function computes the jaccard similarity between the
           % models based on the reaction presence. 
           
           
           integrated_slot_matrix = feval(slot, obj);   % call the function with arguments

           J = squareform(pdist(integrated_slot_matrix','jaccard'));
           disp("Jaccard similarity:")
           1-J
           fig = plot_clustergram(1-J,...
                     obj.model_names,...
                     obj.model_names,...
                     {'Model similarity based on Jaccard distance of rxns existence in the model!'},...
                     [100 100 800 600]);   
        end
        
        function [fig,J] = get_jaccard_similarity_ess_genes(obj,threshold)
           % this function computes the jaccard similarity between the
           % models based on the reaction presence. 
            

           J = squareform(pdist(obj.gene_essentiality.ratio' > threshold,'jaccard'));
           disp("Jaccard similarity:")
           1-J
           fig = plot_clustergram(1-J,...
                     obj.model_names,...
                     obj.model_names,...
                     {'Model similarity based on Jaccard distance of rxns existence in the model!'},...
                     [100 100 800 600]);   
        end
        
        function idx = get_intersection_plot(obj,slot_name,models_to_compare)
            % this function returns the intersection size between all the
            % models + a venn diagramm visualizing the intersection size
            % intersection of up to 4 models
            arguments
               obj 
               slot_name
               models_to_compare (1,:) string = string(obj.model_names(1:min(4,end)))
            end
            [~,idx] = ismember(string(obj.model_names), models_to_compare);
            models_to_compare = models_to_compare(idx(idx ~= 0));
            
            if numel(models_to_compare) > 4
                error('models_to_compare can have at most 4 elements.');
            end
            
            M = feval(slot_name, obj);  
            [~, idx] = ismember(models_to_compare,string(obj.model_names));
            M = M(:,idx);
            idx = plot_flexible_venn(M, models_to_compare);
        end
        function [obj] = get_pathway_prescence(obj,exp,slot)
            arguments
                    obj
                    exp
                    slot ="reaction_presence"
            end
            
           
            % get pathway ids in model 
            set_labels = string(obj.model_names)';
            M = feval(slot, obj); 
            pathways = string(exp.original_model.subSystems);
            unique_pathways = unique(pathways);

            % For each unique pathway, find row indices where it occurs
            groups = arrayfun(@(x) find(pathways == x), unique_pathways, 'UniformOutput', false);
            num_rows = size(M,1);
            num_groups = length(groups);
            G = zeros(num_groups, num_rows);

            for g = 1:num_groups
                G(g, groups{g}) = 1;
            end

            group_sums = G * M;
            % Create table
            T = array2table(group_sums, 'VariableNames', set_labels);

            % Assign row names (only works if unique_pathways is a cellstr or string array)
            T.Properties.RowNames = cellstr(unique_pathways); 
            
            % groupcounts & unique both sort the same way/ therefore can be
            % assigned without matching
            T.original = groupcounts(pathways);
            
            obj.pathway_counts = T
        end

        function scatter_plot_pathway_presence(obj,models_to_compare,plotting_data,labeling_data)

        n = numel(models_to_compare);

        figure
        for i = 1:n
            for j = 1:n
                if i <= j
                    % skip diagonal and upper triangle
                    continue
                end
                
                % create subplot at the correct position
                subplot(n-1, n-1, ((j)+(i-2)*(n-1)))
                
                x = plotting_data.(models_to_compare(j));
                y = plotting_data.(models_to_compare(i));
                scatter(x, y, 20, 'filled')
                axis tight
                
                if i == n
                    xlabel(models_to_compare(j))
                end
                if j == 1
                    ylabel(models_to_compare(i))
                end
                % add text labels for selected pathways
                text(labeling_data.(models_to_compare(j)), labeling_data.(models_to_compare(i)), ...
                     labeling_data.Properties.RowNames, ...
                     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
                     'FontSize', 8, 'FontWeight', 'bold')
            end
        end
        end

        function joined_output = join_fba_output(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            %all_rxns = cellfun(@(x) obj.condition_models.(x).rxns, string(fieldnames(obj.condition_models)),'UniformOutput',false);
            %all_rxns = unique(vertcat(all_rxns{:}));
            all_rxns = obj.fastcore_models.original_model.rxns;

            samples_ordered = arrayfun(@(x) get_sampling_orig_order(obj.fastcore_models.condition_models.(x),obj.fba.(x).v,all_rxns), ...
                                              string(fieldnames(obj.fastcore_models.condition_models)),...
                                              'UniformOutput',false);

            
            biomass_idx = find(ismember(all_rxns, "biomass_reaction"))

            fba_solution = array2table(cell2mat(samples_ordered'), 'VariableNames', fieldnames(obj.fastcore_models.condition_models));
            fba_solution.Properties.RowNames = all_rxns;
            %fba_solution = fba_solution(obj.original_model.rxns,:);
            joined_output = fba_solution;
        end

        function [obj,fluxsum] = compute_flux_sum(obj,slot, figflag,compute_based_on_incoming_flux)
            %COMPUTE_FLUX_SUM this function calculates the fluxsum based on all the
            %rxns producing a metabolite, using the sampling data and the stochiometric
            %matrix from the model!
            
            arguments
               obj (1,1) 
               slot (1,1) ="join_fba_output"
               figflag (1,1) double {mustBeMember(figflag,[1,0])} =1
               compute_based_on_incoming_flux (1,1) double {mustBeMember(compute_based_on_incoming_flux,[1,0])} =1
            end

            solutions = feval(slot, obj); 
            stochiomet = obj.fastcore_models.original_model.S;
            fluxsum=zeros(size(stochiomet,1),size(solutions,2));
            for counter=1:size(solutions,2)
                v=solutions{:,counter}; % one sample
                temp=repmat(v',size(stochiomet,1),1); %
                fluxes=stochiomet.*temp;
                if compute_based_on_incoming_flux
                    fluxSumP=full(sum((fluxes>0).*fluxes,2));
                else
                    fluxSumP=full(sum((fluxes<0).*fluxes,2));
                end
                
                fluxsum(:,counter)=fluxSumP;
            end
            disp('... fluxSum calculated ...')
            
            
            [~,b] = sort(var(fluxsum,0,2),'descend'); % compute the variance over the samples per metabolite
            
            top30_fluxsum = fluxsum(b(1:30),:);
            top30_met_names = obj.fastcore_models.original_model.mets(b(1:30));
            
            if figflag
                figure
                h = boxplot(top30_fluxsum','Labels',top30_met_names);
                set(gca,'FontSize',10,'XTickLabelRotation',45)
                title("30 metabolites with the highest variance for the fluxsum in the samples")
                hold on
            
                numMet = size(top30_fluxsum,1);
                numSamples = size(top30_fluxsum,2);
            
                % Define colors for samples
                colors = lines(numSamples); 
            
                % Sample names for legend
                sampleNames = arrayfun(@(x) sprintf('Sample %d',x), 1:numSamples, 'UniformOutput', false);
            
                % Get x positions of boxes
                bx = findobj(gca,'Tag','Box'); % get all boxes
                xCenters = zeros(numMet,1);
                for k = 1:numMet
                    verts = get(bx(numMet-k+1),'XData'); 
                    xCenters(k) = mean(verts(1:4)); % center of each box
                end
            
                % Overlay scatter points and store handles for legend
                scatterHandles = gobjects(numSamples,1); % initialize
                for j = 1:numSamples
                    for i = 1:numMet
                        scatterHandles(j) = scatter(xCenters(i) + (j-1)*0.05 - 0.025, top30_fluxsum(i,j), ...
                                                    50, colors(j,:), 'filled', 'MarkerEdgeColor','k');
                    end
                end
            
                % Add legend
                legend(scatterHandles, string(obj.model_names), 'Location','best')
                hold off
            end

            
        end


        function get_essentiallity_plots(obj,file_path_figure)
            measures = ["rate", "ratio"];  % list of measures to loop over

            for i = 1:length(measures)
                measure = measures(i);
                %measure = "rate" % "ratio"
                data = obj.gene_essentiality.(measure);

                if measure == "ratio"
                    ylabel_text = 'growth rate KO/WT';
                elseif measure == "rate"
                    ylabel_text = 'growth rate KO';
                end

                % Sort each column individually
                sorted_data = sort(data, 1);  % sorts each column independently
                tol = 1e-6;  % adjust if needed
                max_vals = max(sorted_data, [], 1);
                rows_to_remove = any(abs(sorted_data - max_vals) < tol, 2);
                sorted_data = sorted_data(~rows_to_remove, :);
                % Plot all columns as separate lines
                figure;
                plot(sorted_data, 'LineWidth', 1.5);
                ylabel(ylabel_text);
                xlabel('genes sorted ascending');
                legend(strrep(obj.model_names, '_', '-'), 'Location', 'northwest');
                saveas(gcf,file_path_figure + "\ess_genes_" + measure + ".png");
                grid on;
                hold off;
            end
            
        end
        function top_data = get_pathway_plot(obj,top_pathways,data_type)
                arguments
                    obj
                    top_pathways =1:20
                    data_type ="relative"
                end
            if data_type == "relative"
                relative_counts = array2table(obj.pathway_counts{:,1:end-1} ./ obj.pathway_counts.original);
                relative_counts.Properties.RowNames = obj.pathway_counts.Properties.RowNames;
                relative_counts.Properties.VariableNames = obj.pathway_counts.Properties.VariableNames(1:end-1);
                data = relative_counts
                value_label = "relative counts of subsystem occurence/original model"
            else
                data = obj.pathway_counts(:,1:end-1)
                data{:,:} = data{:,:}./1000
                value_label = "absolute count of rxns per subsystems [# 1/1000]"
            end
            
            
            % Compute variance along rows (dim=2)
            row_var = var(data{:,:}, 0, 2);
            % Get indices of top n highest variance rows
            [~, sortedIdx] = sort(row_var, 'descend');
            
            
            if size(data,1) < length(top_pathways)
                top_pathways = 1:size(data,1)
            end
            top20Idx = sortedIdx(top_pathways);
            
            fig = plot_clustergram(data{top20Idx,:},...
                     data.Properties.RowNames(top20Idx),...
                     data.Properties.VariableNames,...
                     {'Model similarity based on Jaccard distance of rxns existence in the model!'},...
                     [100 100 800 600],...
                     value_label);   
                 
            top_data = data(top20Idx,:);
        end
        function J = get_fba_plot(exp,fba_flux_matrix)
            J = squareform(pdist(fba_flux_matrix{:,:}','jaccard'));
            altcolor =[255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; 
            fig = plot_clustergram(log(1-J),...
                                 string(exp.model_names)',...
                                 string(exp.model_names)',...
                                 {'Similarity of optimal fluxes obtained via FBA [log(jaccard  similarity score)]'},...
                                 [100 100 800 600],...
                                 altcolor);
            
        end
        function J = get_fva_sim_plot(exp,fva_sim_matrix)
            altcolor =[255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; 
            fig = plot_clustergram(fva_sim_matrix,...
                                 string(exp.model_names)',...
                                 string(exp.model_names)',...
                                 {'Similarity of optimal fluxes obtained via FBA [log(jaccard  similarity score)]'},...
                                 [100 100 800 600],...
                                 altcolor);
            disp(fva_sim_matrix)
            
        end
        function obj = add_essentiality_analysis_to_analysis_obj(obj,exp,ratio,rate,rxn_genes,grRateWT)
            
            all_genes = exp.original_model.genes;
           
            fluxsum = arrayfun(@(x) get_order_from_orig_model(exp.condition_models.(x).genes,...
                                                                         ratio.(x),...
                                                                         all_genes,1),...
                                       string(obj.model_names),...
                                       'UniformOutput',false);
                                   
            ess.ratio = cell2mat(fluxsum');
            
            fluxsum2 = arrayfun(@(x) get_order_from_orig_model(exp.condition_models.(x).genes,...
                                                                         rate.(x),...
                                                                         all_genes,grRateWT.(x)),...
                                       string(obj.model_names),...
                                       'UniformOutput',false);
                                   
            ess.rate = cell2mat(fluxsum2');
            
            fluxsum3 = arrayfun(@(x) get_order_from_orig_model(exp.condition_models.(x).genes,...
                                                                         rxn_genes.(x),...
                                                                         all_genes),...
                                       string(obj.model_names),...
                                       'UniformOutput',false);
                                   
            ess.del_rxns = [fluxsum3{:,:}];
            
            obj.gene_essentiality = ess;
            obj.gene_essentiality.grRateWT = grRateWT;
            
        end
        
        function visualize_enrichment(obj)
            
            enrichment = cell2mat(struct2array(structfun(@(x) x.enrichment,obj.enrichment,'UniformOutput',false)));

            choose_gene_sets = find((abs(min(enrichment') - max(enrichment')) > 0.4)' | (sum(enrichment,2)>0.6));
            gene_set_names = arrayfun(@(x) regexprep(x,"_","\_"),obj.enrichment.(obj.model_names{1}).("Database/website"));

            plot_clustergram(enrichment(choose_gene_sets,:),...
                                 gene_set_names(choose_gene_sets)',string(obj.model_names),...
                                 {'enrichment of essential genes per model in gene sets'},...
                                 [100 100 800 600]);

        end
        
        function [fig] = plot_clustergram(obj,data,rownames, colnames,title,position,colorbarLabel,altcolor)
            arguments
                obj
                data
                rownames
                colnames
                title
                position
                colorbarLabel = "Value"  % Default label
                altcolor =[255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; 
            end
          cgo_J = clustergram(data,...
                      'RowLabels', rownames,...
                      'ColumnLabels', colnames,...
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
           cb.Label.FontSize = 12;  % optional formatting

           fig = gcf;
           fig.Position = position;



        end
    end
end



 function rev_find = get_feature_presence(length,indices)
            rev_find = zeros(1,length)';
            rev_find(indices) = 1; 
 end
 
 function [fig] = plot_clustergram(data,rownames, colnames,title,position,colorbarLabel,altcolor)
            arguments
                data
                rownames
                colnames
                title
                position
                colorbarLabel = "Value"  % Default label
                altcolor =[255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; 
            end
          cgo_J = clustergram(data,...
                      'RowLabels', rownames,...
                      'ColumnLabels', colnames,...
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
           cb.Label.FontSize = 12;  % optional formatting

           fig = gcf;
           fig.Position = position;



 end
 
 function idx = plot_flexible_venn(M, set_names)
% M: binary matrix (rows = items, cols = sets)
% set_names: cell array of strings (e.g., {'A','B','C','D'})

n = size(M, 2);  % Number of sets
if n < 2 || n > 4
    error('Only supports 2 to 4 sets.');
end

% Generate all exclusive intersection combinations
labels = strings(2^n - 1, 1);  % Max 15 for 4 sets
idx = struct();  % Max 15 for 4 sets

% Loop through 1 to 2^n - 1 to get all combinations (except all zeros)
for i = 1:(2^n - 1)
    mask = dec2bin(i, n) == '1';  % Logical mask for active sets
    name_set = strjoin(set_names(find(mask)), "_")
    rows_match = all(M(:, mask) == 1, 2) & all(M(:, ~mask) == 0, 2);
    idx.(name_set) = find(rows_match);
    labels(i) = string(sum(rows_match));
end

% Pick colors
% default_colors = [
%     1 0 0;  % Red
%     0 1 0;  % Green
%     0 0 1;  % Blue
%     1 1 0   % Yellow
% ];
% colors = default_colors(1:n, :);

% Call your venn() function
venn(n, ...
    'sets', set_names, ...
    'labels', labels, ...
    'alpha', 0.5, ...
    'edgeC', [1 1 1], ...
    'edgeW', 2);

 end


 function vennfig = venn(n,varargin)
%% Draw venn diagram with two to four sets with optional text labels.
% User can specify the number of sets to draw (maximum four) and label each
% set and the intersectional regions between sets.
% Man Ho Wong (2022).
%
% Input : n [positive integer]
%           Number of sets to draw
%         sets [string | char | cellstr | numeric]
%              An array of set names in left-to-right order
%         labels [string | char | cellstr | numeric]
%                An array of label names for labeling each section;
%                Elements in the array must follow the following order: 
%                For diagram with Set A and B, labels for 3 sections are
%                A, B and A&B.
%                For diagram with Set A, B and C, labels for 7 sections are
%                A, B, C, D, A&B, A&C, B&C and A&B&C.                 
%                For diagram with Set A, B, C and D, labels for 15 sections
%                are A, B, C, D, A&B, A&C, A&D, B&C, B&D, C&D, A&B&C, A&B&D 
%                , A&C&D, B&C&D, A&B&C&D.
%                Any extra labels will be ignored.
%         colors [rows of RGB triplet]
%                Color map for fill colors in left-to-right order.
%                e.g. [1 0 0; 0 1 0; 0 0 1] represents red, green, blue;
%                If number of colors is less than n, colors will be
%                repeated.
%         alpha [0 to 1]
%               Fill color transparency; 0 = fully transparent.
%         edgeC [RGB triplet]
%               Edge color (only effective when 'edgeW' is > 0).
%         edgeW [positive number]
%               Edge width (By default, there is no edge)
%         labelC [RGB triplet]
%                Color of section labels.
%
% Output : A Veenn diagram will be drawn on a new figure.
%          vennfig (optional): A handle to the figure.
%
% Examples: see README.md

% default set names
s = repmat(" ",4,1);   % white space as spaceholder
% default labels
v = repmat(" ",15,1);  % white space as spaceholder
cmap = lines(4);  % color map

% validation functions

%validColor = @(x) ~isempty(validatecolor(x));
%validColors = @(x) ~isempty(validatecolor(x,'multiple')) || validColor;

validColor = @(x) ischar(x) || isstring(x) || ...
    (isnumeric(x) && isvector(x) && length(x) == 3 && all(x >= 0) && all(x <= 1));

validColors = @(x) ischar(x) || isstring(x) || ...
    (isnumeric(x) && size(x,2) == 3 && all(x(:) >= 0) && all(x(:) <= 1));
validNum = @(x) isnumeric(x) && isscalar(x);
validPosNum = @(x) validNum(x) && (x>0);
validPosFrc = @(x) validNum(x) && (x>=0) && (x<=1);

% Build input parser
p = inputParser;
addParameter(p,'sets',s);
addParameter(p,'labels',v);
addParameter(p,'colors',cmap,validColors);
addParameter(p,'alpha', 0.3, validPosFrc);
addParameter(p,'edgeC', 'w', validColor);
addParameter(p,'edgeW', [], validPosNum);
addParameter(p,'labelC', 'k', validColor);


% Parse input
parse(p,varargin{:});
sets = p.Results.sets;
labels = p.Results.labels;
colors = p.Results.colors;
alpha = p.Results.alpha;
edgeC = p.Results.edgeC;
edgeW = p.Results.edgeW;
labelC = p.Results.labelC;

% repeat colors if number of colors given is less than n
%if height(colors) < n
 %   colors = repmat(colors,n/height(colors),1);
%end
numColors = size(colors, 1);
if numColors < n
    colors = repmat(colors, ceil(n / numColors), 1);
end
colors = colors(1:n, :);  % Ensure exactly n rows

% replace spaceholders in f and v with user inputs
%   if user didn't provide enough labels, spaceholders will remain
%   if user provided more labels than needed, extra labels will be ignored
%   thus, the function accepts any number of inputs without causing errors
fRange = min([4 length(sets)]);
for i = 1:fRange
    s(i) = string(sets(i));
end
vRange = min([15 length(labels)]);
for i = 1:vRange
    v(i) = string(labels(i));
end

% for code readability, assign v to variables named by letters
switch n
    case 2
        A = v(1);
        B = v(2);
        AB = v(3);
    case 3
        A = v(1);
        B = v(2);
        C = v(3);
        AB = v(4);
        AC = v(5);
        BC = v(6);
        ABC = v(7);
    case 4
        A = v(1);
        B = v(2);
        C = v(3);
        D = v(4);
        
        AB = v(5);
        AC = v(6);
        AD = v(7);
        BC = v(8);
        BD = v(9);
        CD = v(10);
        
        ABC = v(11);
        ABD = v(12);
        ACD = v(13);
        BCD = v(14);
        
        ABCD = v(15);
end

% figure settings
vennfig = figure('Position',[20 20 800 450],'Color','w');
axis off
daspect([1,1,1])


% circle location and radius
X=1;
Y=1;
r=1;

% draw venn diagram based on number of sets
switch n
    case 2
        xlim([-0.5 4])
        circle(X,Y,r,colors(1,:),alpha);
        circle(X+r,Y,r,colors(2,:),alpha);
        % draw circle A edge again (so it's not covered by circle B)
        circle(X,Y,r,[0 0 0],0);        

        text(1,2.2,s(1),'HorizontalAlignment','right');
        text(2,2.2,s(2),'HorizontalAlignment','left');

        text(0.5,1,A,'HorizontalAlignment','center')
        text(2.5,1,B,'HorizontalAlignment','center')
        text(1.5,1,AB,'HorizontalAlignment','center')

    case 3
        xlim([-0.5 4])
        circle(X,Y,r,colors(1,:),alpha);
        circle(X+r,Y,r,colors(2,:),alpha);
        circle(X+r/2,Y+r,r,colors(3,:),alpha);
        % draw circle A and B edge again (so they are not covered by circle C)
        circle(X,Y,r,colors(1,:),alpha);
        circle(X+r,Y,r,colors(2,:),alpha);

        text(1.5,3.2,s(1),'HorizontalAlignment','center')
        text(-0.1,1,s(2),'HorizontalAlignment','right')
        text(3.1,1,s(3),'HorizontalAlignment','left')

        text(1.5,2.4,A,'HorizontalAlignment','center')
        text(0.5,1,B,'HorizontalAlignment','center')
        text(2.5,1,C,'HorizontalAlignment','center')

        text(1,1.75,AB,'HorizontalAlignment','center')
        text(1.5,0.75,BC,'HorizontalAlignment','center')
        text(2,1.75,AC,'HorizontalAlignment','center')
        
        text(1.5,1.4,ABC,'HorizontalAlignment','center')

    case 4        
        xlim([-3.5 4])

        % ellipse A and B
        [X,Y] = getEllipse(0.8,1.6,[-1.1 1]);
        patch(X,Y,colors(1,:),'FaceAlpha',alpha,'LineStyle','none');
        patch(X+1,Y+0.5,colors(2,:),'FaceAlpha',alpha,'LineStyle','none');

        % ellipse C and D
        [X,Y] = getEllipse(1.6,0.8,[1.1 1]);
        patch(X-1,Y+0.5,colors(3,:),'FaceAlpha',alpha,'LineStyle','none');
        patch(X,Y,colors(4,:),'FaceAlpha',alpha,'LineStyle','none');
        
        % draw ellipse edges separately (so they are not covered by others)
        patch(X-1,Y+0.5,'w','FaceAlpha',0,'LineStyle','none');  % ellipse C
        patch(X,Y,'w','FaceAlpha',0,'LineStyle','none');  % ellipse D
        [X,Y] = getEllipse(0.8,1.6,[-1.1 1]);
        patch(X,Y,'w','FaceAlpha',0,'LineStyle','none');  % ellipse A
        patch(X+1,Y+0.5,'w','FaceAlpha',0,'LineStyle','none');  % ellipse B

        text(-3,3,s(1),'HorizontalAlignment','right')
        text(-2,3.5,s(2),'HorizontalAlignment','right')
        text(2,3.5,s(3),'HorizontalAlignment','left')
        text(3,3,s(4),'HorizontalAlignment','left')
        
        text(-2,1.5,A,'HorizontalAlignment','center')
        text(2,1.5,D,'HorizontalAlignment','center')
        text(-1,2.75,B,'HorizontalAlignment','center')
        text(1,2.75,C,'HorizontalAlignment','center')
        
        
        text(-1.4,2.25,AB,'HorizontalAlignment','center')
        text(1.4,2.25,CD,'HorizontalAlignment','center')
        text(0,2.25,BC,'HorizontalAlignment','center')
        text(-1.25,0.5,AC,'HorizontalAlignment','center')
        text(1.25,0.5,BD,'HorizontalAlignment','center')
        text(0,-0.4,AD,'HorizontalAlignment','center')
        
        text(-0.75,1.5,ABC,'HorizontalAlignment','center')
        text(0.75,1.5,BCD,'HorizontalAlignment','center')
        text(-0.4,0.05,ACD,'HorizontalAlignment','center')
        text(0.4,0.05,ABD,'HorizontalAlignment','center')
        
        text(0,0.5,ABCD,'HorizontalAlignment','center')
        
    otherwise
        disp('n must be an integer between 2 and 4.')
end

% Get all text objects
h=vennfig.findobj('Type','text');

% Configure texts
set(h,'fontsize',11,'FontWeight','bold');
for i = 1:length(h)
    if ismember(h(i).String,sets)
        h(i).FontSize = 14;
        h(i).FontWeight = 'bold';
    else
        h(i).Color = labelC;
    end
end

% Configure edges
if n > 3
    obj = 'patch';
else
    obj = 'rectangle';
end
h=vennfig.findobj('Type',obj);
set(h,'EdgeColor',edgeC);
if ~isempty(edgeW)
    set(h,'LineStyle','-');
    set(h,'LineWidth',edgeW);
end

%%
function [x,y] = getEllipse(r1,r2,C)
beta = linspace(0,2*pi,100);
x = r1*cos(beta) - r2*sin(beta);
y = r1*cos(beta) + r2*sin(beta);
x = x + C(1,1);
y = y + C(1,2);
end

%%
function circle(cX,cY,r,faceC,alpha)
x = cX-r;
y = cY-r;
d = 2*r;
fC = [faceC alpha];
rectangle('Position',[x y d d],'Curvature',1,'FaceColor',fC,'LineStyle','none');
end

 end

function [sampling_fluxsum_ordered] = get_order_from_orig_model(m,s,mets_all,fill_value)
            arguments
               m
               s
               mets_all
               fill_value =0
            end
                                        [~,mapping_mets_in_orig_idx] = ismember(m,mets_all);
                                        if isa(s, 'double')  % or replace with actual numeric class if needed
                                           
                                            sampling_fluxsum_values = ones(length(mets_all), size(s, 2))*fill_value;
                                        elseif isa(s, 'cell')
                                            sampling_fluxsum_values = cell(length(mets_all), size(s, 2));
                                        else
                                            error('Unsupported class for variable s: %s', class(s));
                                        end
                                        sampling_fluxsum_values(mapping_mets_in_orig_idx,:) = s;
                                        sampling_fluxsum_ordered = sampling_fluxsum_values;
                               
end

function [sampling_ordered] = get_sampling_orig_order(m,s,rxns_orig)
                                [~,mapping_rxns_in_orig_idx] = ismember(m.rxns,rxns_orig);
                                sampling_values = zeros(length(rxns_orig),size(s,2));
                                sampling_values(mapping_rxns_in_orig_idx,:) = s;
                                sampling_ordered = sampling_values;
%                                 id_biomass_ordered = find(matches(rxns_orig,"biomass_reaction"));
%                                 id_biomass = find(matches(m.rxns,"biomass_reaction"));
                                
%                                 min(sampling_values(id_biomass_ordered,:))
%                                 max(sampling_values(id_biomass_ordered,:))
%                                 
%                                 min(s(id_biomass,:))
%                                 max(s(id_biomass,:))
                                
end

function [enrichment] = GeneEnrichments(GeneList,GeneList_all)
% gene list should be Symbols

M = numel(GeneList_all);
N = numel(GeneList);

load('cancer_genes_datamined.mat')
K=[];
x=[];
for i=1:numel(cancer_genes_names)
    source = lower(cancer_genes{i});
    
    K = [K;sum(ismember(lower(GeneList_all), source))];
    x = [x;sum(ismember(lower(GeneList), source))];
end


enrichment = table(cancer_genes_names,num2cell(1-hygecdf(x-1,M,K,N)),num2cell(repmat(M,822,1)),num2cell(K),num2cell(repmat(N,822,1)),num2cell(x),...
    'VariableNames',{'Database/website','enrichment','Recon genes','Recon essential genes','Genelist','GeneList essential'});

% [B, I] = sort(1-hygecdf(x,M,K,N),'ascend');
% 
% enrichment = enrichment(I,:);
end





