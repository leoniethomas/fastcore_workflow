## Analysis Pipeline - Documentation

Our analysis pipeline has been divided in three main steps: 

1. Model Building & Project Initialization
2. Single Model Analysis
3. Model Comparison

### 1. Model Building & Project Initialization

The pipeline has been designed for models built with rFastcormics. However, fields specifically associated with rFastcormics being optional, the pipeline can be used for any kind of COBRA-model. 
The pipeline relies on a unique object that we called `project`, following the below format:

#### Model field organization inside project
```text
project
└── models
    └── Name1
        ├── model
        ├── sampleMetadata
        ├── discretizedData
        ├── expressionData
        ├── coreReactions
        └── settings
            ├── medium
            │   ├── mediumComposition
            │   └── manuallySetBoundaries
            │       ├── closedImports
            │       ├── closedExports
            │       ├── unconstrainedImports
            │       └── unconstrainedExports
            ├── scriptParameters
            │   ├── consensusProportion
            │   └── sampleLabeling
            ├── dico
            ├── objFunction
            ├── referenceModel
            ├── mapping
            └── optionalSettings
                ├── medium
                ├── notMediumConstrained
                └── func
```

Project initialization should be done using the `createProject` function. A manual initialization can also be done as long as the above format is respected. 
New models can be added to an already existing project using the `addModelsToProject.m` function. Any new model field will follow the architecture presented above for model `Name1`.

A complete tutorial on how to build models using rFastcormics and initialize a project is available through the `no1_modelBuildingAndProjectInit.m` script. 

### 2. Single Model Analysis

### 3. Model Comparison
