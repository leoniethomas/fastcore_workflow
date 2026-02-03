%% Minimal ACHR sampling null test in COBRA

% Load a model
initCobraToolbox(false)  % Initialize COBRA (no updates)
%model = readCbModel('ecoli_core_model.mat'); % example COBRA model
% load the project object 
working_path = "/Users/leonie.thomas/Documents/fastcore_workflow_with_vanille";
cd (working_path)
addpath(genpath(working_path))

load(working_path + filesep + "context_specific_models" + filesep + "20260119_1042" + filesep + "project_23012026_1453_28012026_1508_obj_vanille_sampling.mat")
model = project.models.WT.model;

%% few things to adjust ? 
% - perform thining, that means throwing out 

%% Set sampling parameters
nSamples = 2000;  % number of samples per chain
nWarmup  = 2000;  % warm-up points
thinStep = 10;    % thinning interval

options.nPointsReturned = 2000;
options.nFiles = 10;
options.maxTime = 36000;
options.nWarmupPoints = 2*size(model.S, 2);
options.nStepsPerPoint = size(model.S, 2)
options.toRound = 1;
options.thinStep = 1000;                % thinning interval (keep every 10th sample)

%% --- First independent ACHR chain ---
rng(1);  % set random seed for reproducibility
sampleFile1 =  char(string(datetime("now", "Format", "yyyyMMdd_HHmm")));
[modelSampling1, samples1] = sampleCbModel(model, sampleFile1,  'ACHR', options, model);

%% --- Second independent ACHR chain ---
rng(2);  % different seed => independent chain
sampleFile2 =  char(string(datetime("now", "Format", "yyyyMMdd_HHmm"))+ ".mat");
[modelSampling2, samples2] = sampleCbModel(model, sampleFile2,  'ACHR', options, model);

%% --- Compare flux distributions ---
% Calculate mean flux per reaction
mean1 = mean(samples1, 1);
mean2 = mean(samples2, 1);

% Example: simple paired t-test per reaction
pValues = zeros(length(mean1), 1);
for i = 1:length(mean1)
    [~, p] = kstest2(samples1(:, i), samples2(:, i));
    pValues(i) = p;
end

%% --- Summary ---
fprintf('Number of "significant" reactions (p < 0.05): %d\n', sum(pValues < 0.05));

% Optional: histogram of p-values
figure;
histogram(pValues, 20);
xlabel('p-value'); ylabel('Number of reactions');
title('ACHR null test: same-model independent chains');

%% do the same with CHRR

%% --- First independent CHRR chain ---
rng(1);  % set random seed for reproducibility
CHRR_sampleFile1 =  char(string(datetime("now", "Format", "yyyyMMdd_HHmm")));
% Preprocess model for CHRR
[CHRR_modelSampling1, CHRR_samples1] = sampleCbModel(model, CHRR_sampleFile1,  'CHRR', options);
%sampleCbModel(model, CHRR_sampleFile1,  'CHRR', options, model);
sample_idx = randsample(1:size(CHRR_samples1,2),2000);
CHRR_samples1 = CHRR_samples1(:,sample_idx);

%% --- Second independent ACHR chain ---
rng(2);  % different seed => independent chain
CHRR_sampleFile2 =  char(string(datetime("now", "Format", "yyyyMMdd_HHmm"))+ ".mat");
[CHRR_modelSampling2, CHRR_samples2] = sampleCbModel(model, CHRR_sampleFile2,  'CHRR', options);
% pick randomly 2000 samples
sample_idx = randsample(1:size(CHRR_samples2,2),2000);
CHRR_samples2 = CHRR_samples2(:,sample_idx);

%% --- Compare flux distributions ---
% Calculate mean flux per reaction
CHRR_mean1 = mean(CHRR_samples1, 1);
CHRR_mean2 = mean(CHRR_samples2, 1);

% Example: simple paired t-test per reaction
CHRR_pValues = zeros(length(CHRR_mean1), 1);
for i = 1:length(CHRR_mean1)
    [~, p] = kstest2(CHRR_samples1(:, i), CHRR_samples2(:, i));
    CHRR_pValues(i) = p;
end



%% --- Summary ---
fprintf('Number of "significant" reactions (p < 0.05): %d\n', sum(pValues < 0.05));

% Optional: histogram of p-values
figure;
histogram(CHRR_pValues, 20);
xlabel('p-value'); ylabel('Number of reactions');
title('CHRR null test: same-model independent chains');


%% Example: two histograms on the same plot
figure;
hold on;

% Plot first histogram
h1 = histogram(CHRR_pValues, 20, 'FaceColor', 'blue', 'FaceAlpha', 0.5);

% Plot second histogram
h2 = histogram(pValues, 20, 'FaceColor', 'red', 'FaceAlpha', 0.5);

xlabel('p-value from kstest2');
ylabel('Number of reactions');
title('P-values comparing two samplings run on the exact same model!');

legend({'chrr', 'achr'});
set(gca, 'FontSize', 20)
hold off;


%% plot the # False positives -> False positive rate

threshold = 0.05;
FPR = sum(pValues < threshold)/size(CHRR_samples2,1);
CHRR_FPR = sum(CHRR_pValues < threshold)/size(CHRR_samples2,1);

bar([1 2], [FPR CHRR_FPR] *100)
ylim([0 100])

xticks([1 2])
xticklabels({'achr','chrr'})
ylabel('False Positive Rate (%)')
title("False Positive Rate in context specific models")
set(gca, 'FontSize', 20)
