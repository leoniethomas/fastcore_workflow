% Loop over each model in parallel fashion
% Working directory
addpath(genpath("/mnt/irisgpfs/users/egonzalez/rFASTCORMICS"))
cd("/mnt/irisgpfs/users/egonzalez/sampling_test_jeff/")
% Loading the models
load('models_medium')
% Configure cluster characteristics
maxWorkers = 8;
cluster = parcluster();
cluster.NumWorkers = maxWorkers;
% Actual parallel loop
parfor i = 1:maxWorkers
    % Make sure each worker initialises CobraToolbox and IBM Cplex solver
    addpath(genpath('/opt/apps/resif/aion/2020a/epyc/software/CPLEX'))
    addpath(genpath("/mnt/irisgpfs/users/egonzalez/cobratoolbox"))
    changeCobraSolver("ibm_cplex")
    pause(30) % wait for 30 secs so ibm cplex can be loaded
    % Set variable to collect results
    samplingResults_medium = struct();
    data = models_medium
    this_model = data(i).models;
    disp(strcat('Hi, I am worker #', num2str(i), ', and I am happy to do your job! Below is my model.'))
    disp(this_model)
    % Set the sampling options common for all models
    options=[];
    options.nPointsReturned = 1000;
    options.nFiles = 10;  % increase this with the nPointsReturned (ratio 1 file ~ 100 samples)
    options.maxTime = 36000;  % 10 hour
    options.nWarmupPoints = 2*size(this_model.S, 2);
    options.nStepsPerPoint = size(this_model.S, 2);
     % Perform sampling
    name=strcat('model_', num2str(i));
    [modelSampling, samples] = sampleCbModel(this_model, name, 'ACHR', options);
    % Store the sampling results in the cell array
    samplingResults_medium(i).samples = samples;
    samplingResults_medium(i).modelSampling = modelSampling;
    parsave_custom(strcat("samplingResults_medium_", num2str(options.nPointsReturned) , "_model_", num2str(i),  ".mat"), samplingResults_medium(i));
end
quit