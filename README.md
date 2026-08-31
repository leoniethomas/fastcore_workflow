# Metabolic network analysis workflow

TackleMMe is based on the creation of one unique MATLAB object, called project. A project is made to store everything, from model building, to model analysis and model comparison.

## 5 steps of the pipeline
  1. **Model building**
  2. **Single-model analysis**
  3. **Models comparison**

Documentation: https://sysbiolux.github.io/analysisPipelineLVT/

## Storing of the data

Everything is stored in a structure named `project`

```
project
├── models
|   ├── orig_model
|   |   └── model
|   ├── consistent_medium_constrained_model
|   |   ├── model
|   |   └── settings
|   |       └── medium
|   ├── context_specific_model_A
|   |   ├── expression_data
|   |   ├── discretized_data
|   |   ├── model
|   |   ├── core_reactions
|   |   ├── settings
|   |   |   ├── reference_model
|   |   |   ├── dico
|   |   |   ├── medium
|   |   |   ├── obj_function
|   |   |   └── optional_settings
|   |   └── analysis_date
|   |       ├── parameters
|   |       ├── FBA
|   |       ├── FVA
|   |       ├── sampling
|   |       ├── single_gene_deletion
|   |       └── double_gene_deletion
|   └── context_specific_model_B
└── comparisons
    └── modelA_vs_modelB
```







