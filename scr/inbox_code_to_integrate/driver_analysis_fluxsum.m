% fluxsum for key metabolites based on FBA for biomass
clearvars -except solverOK, clc, close all

% list models to be analyzed
modelDir=''
models={'mfc_High.mat','mfc_Low.mat'}
names={'High','Low'}

sampleNrs=1:10
modelDir='models\'
models={'model_sample1.mat','model_sample2.mat','model_sample3.mat','model_sample4.mat','model_sample5.mat'}
models = arrayfun(@(x) ['model_sample' num2str(x) '.mat'], sampleNrs, 'UniformOutput', false)
names = arrayfun(@(x) ['S' num2str(x)], sampleNrs, 'UniformOutput', false)

% load generic model
load medium_constrained_consistent_model
model=medium_constrained_consistent_model
% unnest subsystems
if length(model.subSystems{1}) == 1
    disp('unnesting subsystems')
    model.subSystems = vertcat(model.subSystems{:});
end

% pathways for which listed metabolites shall be included
nrPathways=[1] %1,2,4,5,6

mets={};
for counterP=1:numel(nrPathways)
    switch nrPathways(counterP)
        case 1
            pathway='Glycolysis'
            metList={'glc_D[c]','g6p[c]','f6p[c]','fdp[c]','dhap[c]','g3p[c]','13dpg[c]','3pg[c]','2pg[c]','pep[c]','pyr[c]'}
        case 2
            pathway='PPP'
            metList={'ru5p_D[c]','s7p[c]','e4p[c]'}
        case 3
            pathway='TCA_(cytoplasm only)'
            metList={'cit[c]','icit[c]','akg[c]','succ[c]','fum[c]','oaa[c]'}
        case 4
            pathway='Tricarboxylic acid cycle'
            metList={'cit[m]','icit[m]','akg[m]','succoa[m]','succ[m]','fum[m]','mal_L[m]','oaa[m]'}
        case 5
            pathway='Oxidative phosphorylation'
            metList={'nadh[m]','fadh2[m]','focytC[m]','q10h2[m]','atp[m]'}
        case 6
            pathway='Fatty acid oxidation'
            metList={'accoa[c]','accoa[m]','coa[c]','coa[m]'}
        case 7
            pathway='Fatty acid oxidation'
            temp=find(ismember(model.subSystems,pathway));
            metList=findMetsFromRxns(model,model.rxns(temp))';
            temp=find(contains(metList,'coa'));
            metList=metList(temp)
            pathway='Fatty acid oxidation__all__coa'
        case 8
            % Fatty acid oxidation: all carnitin metabolites
            pathway='Fatty acid oxidation'
            temp=find(ismember(model.subSystems,pathway));
            metList=findMetsFromRxns(model,model.rxns(temp))';
            temp=find(contains(metList,'crn'));
            metList=metList(temp)
            pathway='Fatty acid oxidation__all__crn'
        case 9
            pathway='Cholesterol'
            metList={'sql[r]','zymst[r]','chsterol[r]','chsterol[c]'}
            %             removed: 'mev_R[c]','5pmev[c]',
        case 10
            %             pathway='Tyrosine metabolism'
            %             temp=find(ismember(consistent_model.subSystems,pathway));
            %             metList=findMetsFromRxns(consistent_model,consistent_model.rxns(temp))';
            pathway='Dopamine / Tyrosine metabolism'
            metList={'phe_L[c]','tyr_L[c]','34dhphe[c]','dopa[c]','34dhpac[c]','34dhpha[c]','34dhpe[c]','tym[c]'}
            % removed:
    end
    mets=[mets,metList];
end
mets

%%
resBio=[];
clear Tres
for counter=1:numel(models)
    counter
    % load model
    load([modelDir, models{counter}])
    m=mfc;
    
    % incease O2 uptake
    m.lb(find(ismember(m.rxns,'EX_o2[e]')))=-10000; %%%%%%%%%%%%Oxygen uptake
    %     m.lb(find(m.lb==-1000))=-10000;
    %     m.ub(find(m.ub==1000))=10000;
    
    % FBA for biomass
    sol = optimizeCbModel(m,'max','zero')
    sol.f
    resBio=[resBio, sol.f];
    
    % calculate flux sum per metabolite
    temp=repmat(sol.x',size(m.S,1),1);
    fluxes=m.S.*temp;
    fluxSumP=full(sum((fluxes>0).*fluxes,2));
    
    % extract fluxsum for listed metabolites
    temp=find(ismember(m.mets,mets));
    T=table(m.mets(temp),fluxSumP(temp),'VariableNames',{'Metabolite',names{counter}});
    if exist('Tres')
        Tres=outerjoin(Tres,T,'Keys','Metabolite','MergeKeys',true);
        temp=table2array(Tres(:,2:end));
        temp(isnan(temp))=0;
        Tres(:,2:end)=array2table(temp);
    else
        Tres=T;
    end
end
% output results
io=array2table(resBio,'VariableNames',names)
Tres=sortrows(Tres,2,'descend')
