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

A complete tutorial on how to run analyses on a single model is available through the `no2_singleModelAnalysis.m` script.

The `singleModelAnalysis` function allows to perform one or several analyses (FBA, FVA, sampling, loopless sampling, KL divergence, single and double gene deletion) on one or multiple models within a project. Results and parameters are stored for each analysis run, identified by a unique timestamp-based ID. A checkpoint system is available to resume analyses in case of a crash.

Individual analyses can be added to an existing run using the `addAnalysisToExistingOne` function, with interactive confirmation when overwriting existing results.

A PDF report summarizing the main results (model characteristics, exchange fluxes, pathway-level fluxes) can be generated using the `writeAnalysisReport` function.

Details about the available analyses, parameter tables, and report generation are available in https://sysbiolux.github.io/analysisPipelineLVT/single_model_analysis/.

### 3. Model Comparison

A complete tutorial on how to compare models is available through the `no3_modelsComparison.m` script.

Before running a comparison, the `chooseActiveAnalysis` function must be called to designate which analysis run to use for each model. The `modelsComparison` function then performs three types of comparison: structural (presence/absence of reactions, metabolites, and genes), functional (FBA fluxes, FVA similarity, pathway enrichment), and sampling-based (flux distributions, inter-model KL divergence). All results and generated plots are stored in a dedicated `comparisons` field within the project.

Details about the comparison types and their outputs are available in https://sysbiolux.github.io/analysisPipelineLVT/model_comparison/.

### Project overview

A complete `project`, as well as intermediate projects and workspaces, can be downloaded here : https://zenodo.org/records/22209352?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjBkN2Y0ZTI3LTU0ODItNDQxYS05NGJhLTY1MWI2NTc4ODhlMCIsImRhdGEiOnt9LCJyYW5kb20iOiI0OWJjMzFmOWY4MGZkOGU4ZjVkZmI4NDY5NDhiZTQ3ZiJ9.ii7lRdgZ1I5C2fUK3J19wJUdFeSRUkJf7Vx-ttHiIR31ihqqXqc0EVUgcjaE7PN3Sw-BdNVWk6VxHlv-WzWm1w.
You can also navigate through a `project` by visiting the [project layout](https://sysbiolux.github.io/analysisPipelineLVT/project_layout/) page of our documentation.
