# Tackling Metabolic Models Exploration

TackleMMe is a pipeline developed in MATLAB and developed for tackling Metabolic Models exploration. It is based on the creation of a unique MATLAB object, called `project`. A `project` is made to store everything, from model building, to model analysis and model comparison.

## 3 steps of the pipeline
  1. **Model building and Project Initialization**
  2. **Single model analysis**
  3. **Model comparison**

A documentation, including installation requirements, project layout description and functions documentation is available here: https://sysbiolux.github.io/analysisPipelineLVT/.

## Storage of the data

Everything is stored in a structure named `project`. The architecture looks like the below tree.
A complete `project` can be downloaded here : https://zenodo.org/records/22209352?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjIyYmUwNDIxLWM2ZWQtNDFlYi1iYmUwLTQ1Y2Q1ZWIxMGVkNyIsImRhdGEiOnt9LCJyYW5kb20iOiIwNDBiMDRiNTM2MzNjODdkYzdjODk0MjQ4ODQyNWM5NiJ9.AaGq3zPqrzxOAdq4AtCyCkYjNYpDWlzp_JkbtYQgpUU5hq-cHXDmPp_BqpGvtHbRVlkrOpXiiNAFa5dTaEjKYQ. The object shows how a project looks like after running the entire pipeline. This one specifically corresponds to the tutorial example on Breast Cancer data.

## Running an example

A running example is available in the BRCAexample folder. Associated data are provided in the data folder.





