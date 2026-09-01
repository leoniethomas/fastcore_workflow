%% Main script n°1: Model Building & Project Initialization
% This script is showing how to build context-specific models using rFastcormics and how to integrate them in a project format suitable 
% for the rest of the analysis pipeline (singleModelAnalysis + modelsComparison).
% rFastcormics scripts are available from the CobraToolbox or from the
% following github repository: https://github.com/sysbiolux/rFASTCORMICS/tree/master/rFASTCORMICS%20for%20RNA-seq%20data/rFASTCORMICS_v2
% BRCA dataset used in this example was retrieved from TCGA using loadTCGAdata.R
% RPMI medium has been downloaded from https://github.com/sysbiolux/MetabolicMediaLibrary/tree/main/media

%% INITIALIZING THE ENVIRONNEMENT
initCobraToolbox();
changeCobraSolver('gurobi');
feature astheightlimit 2000;

%% LOADING EXAMPLE WORKSPACE
brca_matrix = readtable("data/brca_matrix.csv", 'ReadRowNames', true, 'VariableNamingRule', 'preserve');
samples_metadata = readtable("data/samples_metadata.csv", 'ReadRowNames', true, 'VariableNamingRule', 'preserve');
gene_metadata = readtable("data/gene_metadata.csv", 'ReadRowNames', true, 'VariableNamingRule', 'preserve');
origModel = load('data/Recon3D.mat', 'Recon3D').Recon3D;
dico = load('data/dico.mat').dico;
medium = readtable('data/RPMI1640.tsv', 'FileType', 'text', 'Delimiter', '\t');

%% RENAMING DICO COLUMNS 
dico.Properties.VariableNames{1} = 'geneIdsInModel';
dico.Properties.VariableNames{2} = 'geneIdsInData';

% Bonus: extra column with geneNames (for pipeline figures)
dico.geneNames = dico.geneIdsInData;

%% DISCRETIZATION
brcaMatrixArray = table2array(brca_matrix); % transform table to array
sampleNames = samples_metadata.barcode'; % sample names
discretized = discretizeFPKM(brcaMatrixArray, sampleNames);

%% GET DISCRETIZED MATRICES PER CONDITION
% Separate control vs tumor
isNormal = strcmp(samples_metadata.sample_type, 'Solid Tissue Normal');
isTumor  = strcmp(samples_metadata.sample_type, 'Primary Tumor');

fprintf('Normal Tissue : %d\n', sum(isNormal));
fprintf('Tumors : %d\n', sum(isTumor));

% Get discretized table for the normal condition
discretizedNormal = discretized(:, isNormal);

%groupcounts(samples_metadata, "ajcc_pathologic_stage")

% Get discretized table per cancer stage
isStageI  = (strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage I') | strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage IA') | strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage IB'));
isStageII = (strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage II') | strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage IIA') | strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage IIB'));
isStageIII = (strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage IIIA') | strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage IIIB') | strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage IIIC'));
isStageIV = (strcmpi(samples_metadata.ajcc_pathologic_stage, 'Stage IV'));

idxTumorStageI = find(isTumor & isStageI);
idxTumorStageII = find(isTumor & isStageII);
idxTumorStageIII = find(isTumor & isStageIII);
idxTumorStageIV = find(isTumor & isStageIV);

fprintf('Stage I   : %d samples\n', numel(idxTumorStageI));
fprintf('Stage II  : %d samples\n', numel(idxTumorStageII));
fprintf('Stage III : %d samples\n', numel(idxTumorStageIII));
fprintf('Stage IV  : %d samples\n', numel(idxTumorStageIV));

% Get discretized tables per condition
discretizedStageI = discretized(:, idxTumorStageI);
discretizedStageII = discretized(:, idxTumorStageII);
discretizedStageIII = discretized(:, idxTumorStageIII);
discretizedStageIV = discretized(:, idxTumorStageIV);

% Get expression data per condition
expDataControl = brcaMatrixArray(:, isNormal);
expDataStageI = brcaMatrixArray(:, idxTumorStageI);
expDataStageII = brcaMatrixArray(:, idxTumorStageII);
expDataStageIII = brcaMatrixArray(:, idxTumorStageIII);
expDataStageIV = brcaMatrixArray(:, idxTumorStageIV);

% Get sample metadata per condition
samplesMetadataControl = samples_metadata(isNormal, :);
samplesMetadataStageI = samples_metadata(idxTumorStageI, :);
samplesMetadataStageII = samples_metadata(idxTumorStageII, :);
samplesMetadataStageIII = samples_metadata(idxTumorStageIII, :);
samplesMetadataStageIV = samples_metadata(idxTumorStageIV, :);

% Add a column for labeling (necessary for labeling wanted later in this exercize)
samplesMetadataControl.sample_labeling = repmat('control', height(samplesMetadataControl), 1);
samplesMetadataStageI.sample_labeling = repmat('stageI', height(samplesMetadataStageI), 1);
samplesMetadataStageII.sample_labeling = repmat('stageII', height(samplesMetadataStageII), 1);
samplesMetadataStageIII.sample_labeling = repmat('stageIII', height(samplesMetadataStageIII), 1);
samplesMetadataStageIV.sample_labeling = repmat('stageIV', height(samplesMetadataStageIV), 1);

%%  CHECKING CONSISTENCY OF THE ORIGINAL MODEL
consistentRxnsBool = fastcc(origModel, 1e-4, 1);

consistentModel = removeRxns(origModel, origModel.rxns(setdiff(1:numel(origModel.rxns), consistentRxnsBool))); % create a consistent model based on the vector A
fastcc(consistentModel, 1e-4, 1);

%% CONSTRAINING MODEL BOUNDS USING MEDIUM CONCENTRATIONS
medium.Concentration_uM = medium.Concentration_M*1e6;
mediumConstrainedModel = changeRxnBounds(consistentModel, medium.ExRxns_Recon3D, -medium.Concentration_uM, 'l');
% Adding specific constraints
mediumConstrainedModel = changeRxnBounds(mediumConstrainedModel, {'EX_h2o2[e]', 'EX_o2s[e]', 'EX_oh1[e]', 'EX_ppi[e]'}, 0, 'l'); % closing uptakes
mediumConstrainedModel = changeRxnBounds(mediumConstrainedModel, {'sink_band[c]', 'EX_o2s[e]', 'EX_oh1[e]'}, 0, 'u'); % closing exports
% Checking model consistency
consistentRxnsBoolAfterMedium = fastcc(mediumConstrainedModel, 1e-4, 1);
consistentMediumConstrainedModel = removeRxns(mediumConstrainedModel, mediumConstrainedModel.rxns(setdiff(1:numel(mediumConstrainedModel.rxns),...
    consistentRxnsBoolAfterMedium)));

consistentMediumConstrainedModel = changeObjective(consistentMediumConstrainedModel, "biomass_reaction");
%% SETTINGS
rownames = gene_metadata.gene_name; % gene_names, needed for rFastcormics
biomassRxn = 'biomass_reaction';
consensusProportion = 0.75;
epsilon = 1e-4;
optionalSettings = {};
optionalSettings.medium = medium.Mets_Recon3D; % the rest will be closed
optionalSettings.notMediumConstrained = {'EX_o2[e]', 'EX_h2o[e]', 'EX_co2[e]', 'EX_h[e]'}; % stay opened
optionalSettings.func = {'biomass_reaction', 'biomass_maintenance', 'DM_atp_c_'}; % force the inclusion in the model
adaptiveScalingFlag = 0;
fillingMediumFlag = 1;

%% RUNNING rFASTCORMICS
% Extra documentation for rFASTCORMICS can be found here: https://github.com/sysbiolux/rFASTCORMICS/tree/master/rFASTCORMICS%20for%20RNA-seq%20data/rFASTCORMICS_v2

[modelControl, retainedRxnsControl, idxCoreReactionsControl, paramsForPipelineControl] = rFastcormicsPipeline(consistentMediumConstrainedModel, discretizedNormal, rownames, dico, biomassRxn, ...
    consensusProportion, epsilon, optionalSettings, fillingMediumFlag, adaptiveScalingFlag); % Normal Tissue

[modelStageI, retainedRxnsStageI, idxCoreReactionsStageI, paramsForPipelineStageI] = rFastcormicsPipeline(consistentMediumConstrainedModel, discretizedStageI, rownames, dico, biomassRxn, ...
    consensusProportion, epsilon, optionalSettings, fillingMediumFlag, adaptiveScalingFlag); % Stage I

[modelStageII, retainedRxnsStageII, idxCoreReactionsStageII, paramsForPipelineStageII] = rFastcormicsPipeline(consistentMediumConstrainedModel, discretizedStageII, rownames, dico, biomassRxn, ...
    consensusProportion, epsilon, optionalSettings, fillingMediumFlag, adaptiveScalingFlag); % Stage II

[modelStageIII, retainedRxnsStageIII, idxCoreReactionsStageIII, paramsForPipelineStageIII] = rFastcormicsPipeline(consistentMediumConstrainedModel, discretizedStageIII, rownames, dico, biomassRxn, ...
    consensusProportion, epsilon, optionalSettings, fillingMediumFlag, adaptiveScalingFlag); % Stage III

[modelStageIV, retainedRxnsStageIV, idxCoreReactionsStageIV, paramsForPipelineStageIV] = rFastcormicsPipeline(consistentMediumConstrainedModel, discretizedStageIV, rownames, dico, biomassRxn, ...
    consensusProportion, epsilon, optionalSettings, fillingMediumFlag, adaptiveScalingFlag); % Stage IV

save('data/contextSpecificModelsBRCA.mat', 'modelControl', 'retainedRxnsControl', 'idxCoreReactionsControl', 'paramsForPipelineControl', ...
    'modelStageI', 'retainedRxnsStageI', 'idxCoreReactionsStageI', 'paramsForPipelineStageI', ...
    'modelStageII', 'retainedRxnsStageII', 'idxCoreReactionsStageII', 'paramsForPipelineStageII', ...
    'modelStageIII', 'retainedRxnsStageIII', 'idxCoreReactionsStageIII', 'paramsForPipelineStageIII', ...
    'modelStageIV', 'retainedRxnsStageIV', 'idxCoreReactionsStageIV', 'paramsForPipelineStageIV');

%load('data/contextSpecificModelsBRCA.mat');

%% CREATE A STRUCTURE FOR PROJECT CREATION
% Model Name
paramsForPipelineControl.modelName = 'Control';
paramsForPipelineStageI.modelName = 'StageI';
paramsForPipelineStageII.modelName = 'StageII';
paramsForPipelineStageIII.modelName = 'StageIII';
paramsForPipelineStageIV.modelName = 'StageIV';
% Expression Data
paramsForPipelineControl.expressionData = expDataControl;
paramsForPipelineStageI.expressionData = expDataStageI;
paramsForPipelineStageII.expressionData = expDataStageII;
paramsForPipelineStageIII.expressionData = expDataStageIII;
paramsForPipelineStageIV.expressionData = expDataStageIV;
% Medium Composition
[paramsForPipelineControl.mediumComposition, ...
    paramsForPipelineStageI.mediumComposition, ...
    paramsForPipelineStageII.mediumComposition, ...
    paramsForPipelineStageIII.mediumComposition, ...
    paramsForPipelineStageIV.mediumComposition] = deal(medium);
% Manually Set Boundaries
manuallySetBoundaries = struct();
manuallySetBoundaries.closedImports = {'EX_h2o2[e]', 'EX_o2s[e]', 'EX_oh1[e]', 'EX_ppi[e]'};
manuallySetBoundaries.closedExports = {'sink_band[c]', 'EX_o2s[e]', 'EX_oh1[e]'};
manuallySetBoundaries.unconstrainedImports = {'EX_h2o[e]', 'EX_h[e]'};
manuallySetBoundaries.unconstrainedExports = {}; 
[paramsForPipelineControl.manuallySetBoundaries, ...
    paramsForPipelineStageI.manuallySetBoundaries, ...
    paramsForPipelineStageII.manuallySetBoundaries, ...
    paramsForPipelineStageIII.manuallySetBoundaries, ...
    paramsForPipelineStageIV.manuallySetBoundaries] = deal(manuallySetBoundaries);
% Sample Metadata
paramsForPipelineControl.sampleMetadata = samplesMetadataControl;
paramsForPipelineStageI.sampleMetadata = samplesMetadataStageI;
paramsForPipelineStageII.sampleMetadata = samplesMetadataStageII;
paramsForPipelineStageIII.sampleMetadata = samplesMetadataStageIII;
paramsForPipelineStageIV.sampleMetadata = samplesMetadataStageIV;
% Sample Labeling (needs to be a column of sampleMetadata)
[paramsForPipelineControl.sampleLabeling, ...
    paramsForPipelineStageI.sampleLabeling, ...
    paramsForPipelineStageII.sampleLabeling, ...
    paramsForPipelineStageIII.sampleLabeling, ...
    paramsForPipelineStageIV.sampleLabeling] = deal('sample_labeling');
% Reference Model
[paramsForPipelineControl.referenceModel, ...
    paramsForPipelineStageI.referenceModel, ...
    paramsForPipelineStageII.referenceModel, ...
    paramsForPipelineStageIII.referenceModel, ...
    paramsForPipelineStageIV.referenceModel] = deal('consistentMediumConstrainedModel');
% Create a cell with all the params structures
paramsForPipeline = {paramsForPipelineControl, ...
    paramsForPipelineStageI, ...
    paramsForPipelineStageII, ...
    paramsForPipelineStageIII, ...
    paramsForPipelineStageIV};

%% CREATE A PROJECT
BRCAProject = createProject(paramsForPipeline);

%% ADD A MODEL TO AN EXISTING PROJECT
paramsConsistentMediumConstrainedModel = struct();
paramsConsistentMediumConstrainedModel.modelName = "consistentMediumConstrainedModel";
paramsConsistentMediumConstrainedModel.contextSpecificModel = consistentMediumConstrainedModel;
paramsConsistentMediumConstrainedModel.dico = dico;
paramsConsistentMediumConstrainedModel.mediumComposition = medium;
BRCAProject = addModelsToProject(BRCAProject, {paramsConsistentMediumConstrainedModel});

save('data/BRCAProjectNo1.mat', 'BRCAProject')

