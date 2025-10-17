%%
clearvars -except solverOK, clc, close all

% modelDir='allmodels_Soumen151025\'
% load([modelDir 'High group.mat'])
% load([modelDir 'Low group.mat'])
% load('mfc_High.mat')
% load('mfc_Low.mat')
% m=mfc
%
% % some stats
% tabulate(m.ub)
% tabulate(m.lb)
%
% tabulate(m.c)
% m.rxns(find(m.c))

%% growth via FBA
% list models to be analyzed
modelDir=''
models={'mfc_High.mat','mfc_Low.mat'}
names={'High','Low'}

sampleNrs=1:40
modelDir='models\'
models={'model_sample1.mat','model_sample2.mat','model_sample3.mat','model_sample4.mat','model_sample5.mat'}
models = arrayfun(@(x) ['model_sample' num2str(x) '.mat'], sampleNrs, 'UniformOutput', false)
names = arrayfun(@(x) ['S' num2str(x)], sampleNrs, 'UniformOutput', false)

resBio=[];
clear TresUpt TresSec
for counter=1:numel(models)
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
    
    % extract uptake rates
    [exRxns, Ex_orgaInd] = findEX_Rxns(m, m.rxns(find(m.c)));
    temp=sol.x(find(ismember(m.rxns,exRxns)));
    T=table(exRxns(temp<0),temp(temp<0),'VariableNames',{'EX',names{counter}});
    sortrows(T,2)
    if exist('TresUpt')
        TresUpt=outerjoin(TresUpt,T,'Keys','EX','MergeKeys',true);
        temp=table2array(TresUpt(:,2:end));
        temp(isnan(temp))=0;
        TresUpt(:,2:end)=array2table(temp);
    else
        TresUpt=T;
    end
    
    %extract secretion rates
    temp=sol.x(find(ismember(m.rxns,exRxns)));
    T=table(exRxns(temp>0),temp(temp>0),'VariableNames',{'EX',names{counter}});
    sortrows(T,2)
    if exist('TresSec')
        TresSec=outerjoin(TresSec,T,'Keys','EX','MergeKeys',true);
        temp=table2array(TresSec(:,2:end));
        temp(isnan(temp))=0;
        TresSec(:,2:end)=array2table(temp);
    else
        TresSec=T;
    end
    
end
% output results
resBio=array2table(resBio,'VariableNames',names)
TresUpt=sortrows(TresUpt,2)
TresSec=sortrows(TresSec,2,'descend')

%% number of rxns and mets
% list models to be analyzed
sampleNrs=1:40
modelDir='models\'
models={'model_sample1.mat','model_sample2.mat','model_sample3.mat','model_sample4.mat','model_sample5.mat'}
models = arrayfun(@(x) ['model_sample' num2str(x) '.mat'], sampleNrs, 'UniformOutput', false)
names = arrayfun(@(x) ['S' num2str(x)], sampleNrs, 'UniformOutput', false)

resRxns=[];
resMets=[];
for counter=1:numel(models)
    % load model
    load([modelDir, models{counter}])
    m=mfc;
    m.lb(find(ismember(m.rxns,'EX_o2[e]')))=-10000; %%%%%%%%%%%%Oxygen uptake
    
    % extract information
    resRxns=[resRxns, numel(m.rxns)];
    resMets=[resMets, numel(m.mets)];
end

% present results
figure
subplot(1,2,1)
boxplot(resRxns)
title('Number of reactions')
subplot(1,2,2)
boxplot(resMets)
title('Number of metabolites')

figure
plot(resRxns,resMets,'*')
xlabel('Number of reactions')
ylabel('Number of metabolites')

%% pathway presence rates
clearvars -except solverOK, clc, close all

% load generic model
load medium_constrained_consistent_model
model=medium_constrained_consistent_model
% unnest subsystems
if length(model.subSystems{1}) == 1
    disp('unnesting subsystems')
    model.subSystems = vertcat(model.subSystems{:});
end
temp=tabulate(model.subSystems);
Tres=table(temp(:,1),cell2mat(temp(:,2)),'VariableNames',{'pathway','generic'})

% list models to be analyzed
sampleNrs=1:40
modelDir='models\'
models={'model_sample1.mat','model_sample2.mat','model_sample3.mat','model_sample4.mat','model_sample5.mat'}
models = arrayfun(@(x) ['model_sample' num2str(x) '.mat'], sampleNrs, 'UniformOutput', false)
names = arrayfun(@(x) ['S' num2str(x)], sampleNrs, 'UniformOutput', false)

for counter=1:numel(models)
    counter
    %load model
    load([modelDir, models{counter}])
    m=mfc;
    m.lb(find(ismember(m.rxns,'EX_o2[e]')))=-10000; %%%%%%%%%%%%Oxygen uptake
    % unnest subsystems
    if length(m.subSystems{1}) == 1
        disp('unnesting subsystems')
        m.subSystems = vertcat(m.subSystems{:});
    end
    % extract number of active reactions per pathway
    temp=tabulate(m.subSystems);
    Tres2=table(temp(:,1),cell2mat(temp(:,2)),'VariableNames',['pathway',names(counter)]);
    Tres=outerjoin(Tres,Tres2,'Keys','pathway','MergeKeys',true);
    temp=table2array(Tres(:,2:end));
    temp(isnan(temp))=0;
    Tres(:,2:end)=array2table(temp);
end
Tres

% some stats: largest differences (absolute / relative)
data=table2array(Tres(:,3:end));

largestDiff=max(data,[],2)-min(data,[],2);
T=sortrows([Tres(:,1:2), table(largestDiff)],3)

largestDiffRel=largestDiff./table2array(Tres(:,2))
T=sortrows([Tres(:,1:2), table(largestDiffRel)],3)
