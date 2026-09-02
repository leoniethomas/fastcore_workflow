# Welcome to TackleMMe

TackleMMe is a pipeline developed in MATLAB designed for tackling Metabolic Models exploration.

## Concept

TackleMMe is based on the creation of one unique MATLAB object, called *project*. A *project* is made to store everything, from model building, to model analysis and model comparison.
Tutorial scripts can be found in the [BRCA example folder of our github](https://github.com/sysbiolux/analysisPipelineLVT/tree/main/BRCAexample).

A complete `project` object after running the entire pipeline (`20262608_BRCAProject.mat`), as well as intermediate objects after running the several steps of the pipeline are available [here](https://zenodo.org/records/22209352?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjIyYmUwNDIxLWM2ZWQtNDFlYi1iYmUwLTQ1Y2Q1ZWIxMGVkNyIsImRhdGEiOnt9LCJyYW5kb20iOiIwNDBiMDRiNTM2MzNjODdkYzdjODk0MjQ4ODQyNWM5NiJ9.AaGq3zPqrzxOAdq4AtCyCkYjNYpDWlzp_JkbtYQgpUU5hq-cHXDmPp_BqpGvtHbRVlkrOpXiiNAFa5dTaEjKYQ). Workspaces associated with the tutorials are also included.

The pipeline and [tutorials](https://github.com/sysbiolux/analysisPipelineLVT/tree/main/BRCAexample) has been organized in three main steps:

1. **Model building and project initialization**

    The tutorial includes how to:

    - build context-specific models from RNA-seq data using [rFASTCORMICS](https://github.com/sysbiolux/rFASTCORMICS/tree/master/rFASTCORMICS%20for%20RNA-seq%20data/rFASTCORMICS_v2),
    - initialize a project,
    - add a model to an already existing project.

    Although TackleMMe offers the possibility to build context-specific models using rFASTCORMICS, any COBRA-format model can be stored inside a *project*. Mandatory fields required for *project* initialization are detailed in [Project Initialization](project_init.md).

2. **Single Model Analysis**

    Once a project is initialized and contains models, our pipeline allows to:
    
    - perform integrated tests on one or several models,
    - store the results as well as the parameters that have been used for each test,
    - update tests with new parameters,
    - generate a PDF report summarizing the results of the main tests.

    More details about the implemented analysis and the report can be found in [Single Model Analysis](single_model_analysis.md) 

3. **Model Comparison**

    



