## Analysis Pipeline - Documentation

Our analysis pipeline has been divided in three main steps: 

1. Model Building & Project Initialization
2. Single Model Analysis
3. Model Comparison

### 1. Model Building & Project Initialization

A complete tutorial on how to build models using rFastcormics and initialize a project is available through the `no1_modelBuildingAndProjectInit.m` script.

The pipeline has been designed for models built with rFastcormics. However, fields specifically associated with rFastcormics being optional, the pipeline can be used for any kind of COBRA-model. Details about the input data format are available in https://sysbiolux.github.io/analysisPipelineLVT/project_init/.

Project initialization should be done using the `createProject` function. New models can be added to an already existing project using the `addModelsToProject.m` function.

### 2. Single Model Analysis

### 3. Model Comparison
