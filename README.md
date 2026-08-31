# Tackling Metabolic Models Exploration

TackleMMe is a pipeline developed in MATLAB and developed for tackling Metabolic Models exploration. It is based on the creation of a unique MATLAB object, called `project`. A `project` is made to store everything, from model building, to model analysis and model comparison.

## 3 steps of the pipeline
  1. **Model building and Project Initialization**
  2. **Single model analysis**
  3. **Model comparison**

A documentation, including installation requirements, project layout description and functions documentation is available here: https://sysbiolux.github.io/analysisPipelineLVT/.

## Storage of the data

Everything is stored in a structure named `project`. The architecture looks like the below tree.
A complete `project` can be downloaded here : https://zenodo.org/records/22209352?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjFmZDcxYzczLThkZmUtNDEzOC1iM2VhLTI3ZWZkZGE4ZjI2MiIsImRhdGEiOnt9LCJyYW5kb20iOiJkNzVkMjg5ZmJjMWExMTIwN2RkMDIxMThmYzE1NWU1MiJ9.6bOIsZ9qMHhVH76KO6F5b0tNX7AvwxR8-gtQcEn1c5cfNNsV9dC_YzCw5Dtoa3hgZ_IfOc6qJVQr1hTUlXW3vw. The object shows how a project looks like after running the entire pipeline. This one specifically corresponds to the tutorial example on Breast Cancer data.

  ### Project architecture

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

## Running an example

A running example is available in the BRCAexample folder. Associated data are provided in the data folder.





