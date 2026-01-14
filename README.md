# Metabolic network analysis workflow

This repository provides a workflow going from the preprocessing of the RNA-Seq data to the analysis and comparison of the context specific metabolic models created by rFastcormics_v2. 

## 5 steps of the pipeline
  1. RNA preprocessing
  2. Quality control
  **3. Model building**
  **4. Single-model analysis**
  **5. Models comparison**

## Storing of the data

Everything is stored in a structure named `project`

```
project/
├── dico/
├── models/
|   ├── orig_model/
|   ├── consistent_medium_constrained_model/
|   ├── context_specific_model_A/
|   |   ├── discretization/
|   |   ├── model/
|   |   ├── reactions/
|   |   ├── core_reactions/
|   |   ├── settings/
|   |   |   ├── medium/
|   |   |   ├── obj_function/
|   |   |   └── optional_settings/
|   |   └── analysis/
|   |   |   ├── FBA/
|   |   |   ├── FVA/
|   |   |   ├── flux_sum/
|   |   |   ├── sampling/
|   |   |   ├── single_gene_deletion/
|   |   |   └── double_gene_deletion/
|   └── context_specific_model_B/
└── comparisons/
    └── modelA_vs_modelB/
```







