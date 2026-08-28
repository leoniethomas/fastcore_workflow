# Welcome to TeckleMMe

TeckleMMe is a pipeline developed in MATLAB designed for tackling Metabolic Models exploration.

## Concept

TeckleMMe is based on the creation of one unique MATLAB object, called *project*. A *project* is made to store everything, from model building, to model analysis and model comparison.

The pipeline has been organized in three main steps:

1. **Model building and project initialization**

    The tutorial includes how to:

    - build context-specific models from RNA-seq data using rFASTCORMICS,
    - initialize a project,
    - add a model to an already existing project.

    Although TeckleMMe offers the possibility to build context-specific models using rFASTCORMICS, any COBRA-format model can be stored inside a *project*. Mandatory fields required for *project* initialization are detailed in [Project format](project_init.md).

2. **Single Model Analysis**

3. **Model Comparison**

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
