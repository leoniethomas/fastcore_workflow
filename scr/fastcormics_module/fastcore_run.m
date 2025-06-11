classdef fastcore_run
    %fastcore_experiments -> summary of the fastcormics worflow
    %   
    
    properties
        model
        sampling
        fluxsum
    end
    
    methods
        function obj = fastcore_run(model,sampling,run_fluxsum)
            arguments
               model
               sampling
               run_fluxsum (1,1) double {mustBeMember(run_fluxsum,[1,0])} =0
            end
            obj.model = model;
            obj.sampling = sampling;
            if run_fluxsum
                disp("...calculating the fluxsum from the sampling data...")
                [~,obj.fluxsum] = obj.compute_flux_sum(1);
            else
                disp("... in the default setting the fluxsum is not calculated when initializing the a fastcore run object ...")
            end
        end
        
        function [obj,fluxsum] = compute_flux_sum(obj,figflag)
            %COMPUTE_FLUX_SUM this function calculates the fluxsum based on all the
            %rxns producing a metabolite, using the sampling data and the stochiometric
            %matrix from the model!
            
            arguments
               obj (1,1) fastcore_run
               figflag (1,1) double {mustBeMember(figflag,[1,0])} =1
            end

            obj.fluxsum=zeros(size(obj.model.S,1),size(obj.sampling,2));
            for counter=1:size(obj.sampling,2)
                v=obj.sampling(:,counter); % one sample
                temp=repmat(v',size(obj.model.S,1),1); %
                fluxes=obj.model.S.*temp;
                fluxSumP=full(sum((fluxes>0).*fluxes,2));
                obj.fluxsum(:,counter)=fluxSumP;
            end
            disp('... fluxSum calculated ...')
            fluxsum = obj.fluxsum;
            
            [~,b] = sort(var(obj.fluxsum,0,2),'descend'); % compute the variance over the samples per metabolite
            
            top30_fluxsum = obj.fluxsum(b(1:30),:);
            top30_met_names = obj.model.metNames(b(1:30));
            
            if figflag
                figure
                boxplot( top30_fluxsum','Labels',top30_met_names)
                set(gca,'FontSize',10,'XTickLabelRotation',45)
                title("30 metabolites with the highest variance for the fluxsum in the samples")
                
            end
        end
    end
end

