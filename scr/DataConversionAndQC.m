%% loading counts and gene lengths
clc, clearvars -except solverOK
counts = readtable('C:\Users\leonie.thomas\20250225_glynn_bulk_metabolic_model\data\bulkRNAseq\raw_counts_with_genelengths.csv');
% Length should be download from biomart pay attention to the right assembly
% https://www.ensembl.org/biomart/martview/86ee26709c0735e01f47aad88718ba5e
%lengthsOrig = readtable('..\..\Data\lengths.txt'); % I do not have the
%length data!!
%lengthsOrig=lengthsOrig(:,2:3);
%Ensembl= counts.ENSEMBL_ID;

data = table2array(counts(:,3:end)); %log2 scaling with pseudocount 1
disp('Number of reads:')
sum(data)
data = log2(data+1); %log2 scalling with pseudocount 1

% addpath(genpath(pwd));
%cond={'EV', 'EV', 'EV10nME224HR_I',  'EV10nME224HR_I','1uG124hr_I', '1uG124hr_I', 'ERBtransduced','ERBtransduced'};

disp(' ')
disp('... data loaded ...')
disp(' ')

%% PCA for QC: low technical variance?
% with mapcaplot function
mapcaplot(counts{:,3:end}');
% with pca function and manual plotting
% [~,score,~,~,explained,~] = pca(data');
% figure
% hold on
% plot(score(1:2,1),score(1:2,2),'r.', 'MarkerSize',20 )
% plot(score(3:4,1),score(3:4,2),'g.', 'MarkerSize',20)
% plot(score(5:6,1),score(5:6,2),'b.' ,'MarkerSize',20)
% plot(score(7:8,1),score(7:8,2),'k.', 'MarkerSize',20)
% 
% title('PCA')
% xlabel([num2str(1), ' component: ', num2str(explained(1))])
% ylabel([num2str(2), ' component: ', num2str(explained(2))])
% legend({'EV','EV10nME224HR_I','EV10nME224HR_I', 'ERBtransduced' },'location','best')

%% Boxplot (data distributions)
% all data
figure
boxplot(data)

% without zeros
data2=data;
data2(data2==0)=NaN;
figure
boxplot(data2)

%% count to TPM conversion
% identify genes with length information
% [I,IA,IB]=intersect(counts.ENSEMBL_ID,table2array(lengthsOrig(:,1)));
% data=(data(IA,:));

lengths=counts{:,2};
temp=counts{:,3:end};
temp=temp./repmat(lengths,1,size(temp,2)); %length normalization
TPM=1000000*temp./repmat(nansum(temp),size(temp,1),1); %scaling

% The sum of TPM for each sample should be 1e6
sum(TPM,1)

mapcaplot(TPM');

%% count to FPKM conversion
lengths=counts{:,2};
temp=counts{:,3:end};
temp=1000000*temp./repmat(nansum(temp),size(temp,1),1); %scaling
FPKM=temp./repmat(lengths,1,size(temp,2)); %length normalization


sum(FPKM)
mapcaplot(FPKM');
%% Densityplot (as used in rFASTCORMICS later)
figure
hold on
lTPM=TPM;
lTPM(lTPM==0)=NaN; %remove zeros for densityplot
lTPM=log2(lTPM); %log2 scaling
for i=1:size(lTPM,2)
    [probability_estimate,xi] = ksdensity(lTPM(:,i));
    plot(xi,probability_estimate,':k','LineWidth',1);
end
hold off

% histogram of sample 1
figure
hist(lTPM(:,1),100)

TPMt=table(table2array(counts(IA,1:2)),TPM);
% writetable(TPMt, '\TPM.txt')


%% Densityplot (as used in rFASTCORMICS later)
figure
hold on
lTPM=counts{:,3:end};
lTPM(lTPM==0)=NaN; %remove zeros for densityplot
lTPM=log2(lTPM); %log2 scaling
for i=1:size(lTPM,2)
    [probability_estimate,xi] = ksdensity(lTPM(:,i));
    plot(xi,probability_estimate,':k','LineWidth',1);
end
hold off

% histogram of sample 1
figure
hist(lTPM(:,1),100)

TPMt=table(table2array(counts(IA,1:2)),counts{:,3:end});
% writetable(TPMt, '\TPM.txt')



%% Densityplot (as used in rFASTCORMICS later)
figure
hold on
lTPM=FPKM;
lTPM(lTPM==0)=NaN; %remove zeros for densityplot
lTPM=log2(lTPM); %log2 scaling
for i=1:size(lTPM,2)
    [probability_estimate,xi] = ksdensity(lTPM(:,i));
    plot(xi,probability_estimate,':k','LineWidth',1);
end
hold off

% histogram of sample 1
figure
hist(lTPM(:,1),100)

