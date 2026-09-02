# Single Model Analysis

## Performing analysis

The `singleModelAnalysis` function runs one or more analyses on a model or a list of models and stores the results in the project structure.

### Prerequisites

!!! warning "Project format"
    The `project` structure must follow the format described in the [Project Layout](project_layout.md) page. Use `createProject` and `addModelsToProject` to build a valid project before running analyses.

### Available analyses

The following analyses are currently implemented:

| Analysis | Key | Description |
|---|---|---|
| Flux Balance Analysis | `FBA` | Computes the optimal flux distribution maximizing a given objective function |
| Flux Variability Analysis | `FVA` | Determines the range of feasible fluxes for each reaction under given constraints |
| Sampling | `sampling` | Uniform random sampling of the feasible flux space |
| Loopless sampling | `loopless` | Sampling with thermodynamic constraints eliminating internal loops |
| KL divergence | `kld` | Computes Kullback-Leibler divergence between sampling distributions |
| Single gene deletion | `singleGeneDeletion` | Evaluates the effect of individually knocking out each gene |
| Double gene deletion | `doubleGeneDeletion` | Evaluates the effect of pairwise gene knockouts |

!!! note "Upcoming analyses"
    The list of available analyses is actively being expanded. The following are planned for integration:

    - **Enrichment** — pathway enrichment analysis on gene deletion results

    Check back for updates in future releases.

### Parameter table

The parameter table is a three-column table (`Parameter`, `Analysis`, `Value`) that configures the settings for each analysis. Each row associates a parameter name with a specific analysis key and its corresponding value. The `Analysis` column determines which analysis the parameter applies to — use `all` for parameters shared across every analysis (e.g. model reference, objective function), or a specific key such as `FBA`, `sampling`, or `singleGeneDeletion` for analysis-specific settings.

The **first two rows are required** and must use `all` as the analysis value:

| Parameter | Analysis | Value | Description |
|---|---|---|---|
| `modelReference` | `all` | e.g. `Recon3D` | Reference model identifier |
| `objFunction` | `all` | e.g. `biomass_reaction` | Objective function reaction ID |

A default parameter table is provided in the `data/` folder of the repository ([`defaultParametersTable.csv`](https://github.com/sysbiolux/analysisPipelineLVT/blob/main/data/defaultParametersTable.csv)) and can be used as a starting point. You can load it in MATLAB and modify only the rows relevant to your analyses:

```matlab
parameterTable = readtable('data/defaultParametersTable.csv');
```

!!! tip "Filtering"
    When running `singleModelAnalysis` or `addAnalysisToExistingOne`, only the rows matching the requested analyses (and `all`) are used. Rows for other analyses are ignored, so you can keep a single comprehensive table without trimming it manually.

### Function signature

```matlab
project = singleModelAnalysis(project, parameterTable, modelList, analyses, saveCheckpoint, resumeFromCheckpoint)
```

#### Input arguments

| Name | Type | Description | Default |
|---|---|---|---|
| `project` | `struct` | Project structure following the format described in [Project Layout](project_layout.md) | required |
| `parameterTable` | `table` | Parameter table defining analysis settings (at minimum the default one) | required |
| `modelList` | `string array` | Names of the models to analyze | `{}` (all models) |
| `analyses` | `cell array` | List of analyses to perform | `{}` (all available) |
| `saveCheckpoint` | `logical` | Whether to save a checkpoint after each model | `true` |
| `resumeFromCheckpoint` | `logical` | Whether to resume from the last saved checkpoint | `false` |

#### Output

| Name | Type | Description |
|---|---|---|
| `project` | `struct` | The input project with an `analysis` field added to each model |

### Analysis IDs

Each analysis run is assigned a unique identifier based on the current timestamp:

```
analysis_20240815_1430
```

The format is `analysis_yyyyMMdd_HHmm`, where:

- `yyyyMMdd` — date (year, month, day)
- `HHmm` — time (hour, minute)

This ID is used as a field name under `project.models.<modelName>.analysis.<id>` and allows multiple analysis runs to coexist on the same model without overwriting previous results. When comparing models later, you can reference a specific analysis by its ID.

### Parameter storage

For each analysis run, the parameter table used is stored alongside the results:

```matlab
project.models.<modelName>.analysis.<id>.parameters
```

This ensures full traceability — you can always retrieve the exact settings that produced a given set of results, even if the parameter table was modified between runs.

### Loopless sampling

Loopless sampling can be run in two ways:

1. **In the same run as sampling** — if both `sampling` and `loopless` are requested, the sampling results from the current run are automatically used as input for the loopless analysis.

2. **Using a previous sampling run** — if only `loopless` is requested (without `sampling`), you must specify which existing sampling results to use via the `samplingToUse` parameter in the parameter table. This parameter should contain the analysis ID of a previously run sampling analysis, in the same order as the models in `modelList`.

!!! warning "Order matters"
    When using `samplingToUse`, the analysis IDs must be listed in the same order as their corresponding models in `modelList`. Mismatched ordering will produce an error.

### Checkpoint system

The function includes a checkpoint mechanism to protect against crashes during long-running analyses:

- **Automatic save** — after each model completes, the current project state and progress index are saved to `singleModelAnalysis_checkpoint.mat` (if `saveCheckpoint` is `true`).
- **Resume** — set `resumeFromCheckpoint = true` to load the last checkpoint and continue from the next model. The function prints which model it is resuming from.
- **Error handling** — if an error occurs during analysis of a model, the function catches it, prints the error message and stack trace, and returns the partial project. You can then relaunch with `resumeFromCheckpoint = true` to continue.
- **Cleanup** — the checkpoint file is automatically deleted once all models have been successfully analyzed.

### Usage example

```matlab
% Run all available analyses on all models
project = singleModelAnalysis(project, parameterTable);

% Run only FBA and sampling on specific models
project = singleModelAnalysis(project, parameterTable, ...
    ["model1", "model2"], {"FBA", "sampling"});

% Resume after a crash
project = singleModelAnalysis(project, parameterTable, ...
    resumeFromCheckpoint = true);
```

### Adding analyses to an existing run

The `addAnalysisToExistingOne` function adds one or more analyses to an already existing analysis field — identified by its analysis ID — without creating a new timestamped entry. This is useful when you want to supplement a previous run with additional analyses (e.g. adding gene deletions to a run that only contained FBA and FVA).

```matlab
project = addAnalysisToExistingOne(project, parameterTable, modelName, analyses, analysisId)
```

#### Input arguments

| Name | Type | Description | Default |
|---|---|---|---|
| `project` | `struct` | Project structure with an existing analysis field | required |
| `parameterTable` | `table` | Parameter table containing settings for the analyses to add | required |
| `modelName` | `string` | Name of the model | required |
| `analyses` | `string` or `cell array` | Analysis key(s) to add (e.g. `"FBA"`, `{"FBA","sampling"}`) | required |
| `analysisId` | `string` | Existing analysis ID to add to (e.g. `analysis_20240815_1430`) | required |

!!! warning "Overwriting behavior"
    - The **parameter table** stored under the analysis ID is updated: only the rows corresponding to the requested analyses are replaced. Other analyses' parameters are preserved.
    - If an analysis **already exists** with the same key, its results are **overwritten** after user confirmation.
    - If no parameters table is found, a warning is issued and the new parameters are stored from scratch.

#### Confirmation prompts

When an analysis already exists, the function compares the stored parameters with the new ones and prompts for confirmation:

- **Identical parameters** — displays the parameters and asks whether to re-run the analysis.
- **Different parameters** — displays a table showing the differing parameters (existing value vs. new value) and asks whether to overwrite.
- **No existing parameters for that test** — asks whether to overwrite the analysis and store the new parameters.

In all cases, answering `n` aborts the operation with an error.

#### Usage example

```matlab
% Add gene deletion analyses to an existing run
project = addAnalysisToExistingOne(project, parameterTable, ...
    "model1", {"singleGeneDeletion", "doubleGeneDeletion"}, ...
    "analysis_20240815_1430");

% Re-run FVA with updated parameters on an existing run
project = addAnalysisToExistingOne(project, parameterTable, ...
    "model1", "FVA", "analysis_20240815_1430");
```

## Analysis Reports

The `writeAnalysisReport` function generates a PDF report summarizing the results of a single analysis run. It requires that **FBA and FVA** have been performed on the model — these are the minimum required analyses.

### Function signature

```matlab
writeAnalysisReport(project, modelName, analysisId, varargin)
```

#### Input arguments

| Name | Type | Description | Default |
|---|---|---|---|
| `project` | `struct` | Project structure containing analysis results | required |
| `modelName` | `string` | Name of the model to report on | required |
| `analysisId` | `string` | Analysis ID (e.g. `analysis_20240815_1430`) | required |
| `path` | `string` | Output directory for the PDF file | `''` (current folder) |
| `pathwaysOfInterest` | `cell array` | Pathway names to include in the report | `{}` (none) |
| `metsOfInterest` | `cell array` | Metabolite IDs to include in the report | `{}` (none) |

!!! warning "FBA and FVA required"
    The report generation relies on FBA fluxes and FVA ranges. If these analyses were not performed, the function will error. Make sure `FBA` and `FVA` are included in the `analyses` list when running `singleModelAnalysis`.

### Report structure

The generated PDF is organized as follows:

#### 1. Title page

Displays the model name and the analysis date, extracted from the analysis ID and formatted as a readable date (e.g. `15 August 2024, 14:30`).

#### 2. Table of contents

Auto-generated from the chapter headings.

#### 3. Analysis parameters

The full parameter table used for the analysis run, providing complete traceability of the settings.

#### 4. Model characteristics

A summary table including:

| Field | Description |
|---|---|
| Model consistency | Whether the model passes `fastcc` consistency check (`YES`/`NO`) |
| Number of reactions | Total reaction count in the model |
| Consistent reactions | Number of reactions passing the consistency check |
| Uptake reactions | Number of exchange reactions importing metabolites |
| GPR rules | Number of reactions with gene-protein-reaction associations |
| Core reactions retained | Number of core reactions from discretization still present in the model (if applicable) |

#### 5. Medium composition

If a medium is defined in `project.models.<modelName>.settings.medium`, a table of the medium composition is included, showing metabolites, exchange reactions, concentrations (in µM), and fluxes (in mmol/gDW/h). The displayed columns adapt to the model reference (Recon3D or HumanGEM).

#### 6. Table of exchangers

All exchange reactions (`EX_` prefix) with their FBA flux, FVA minimum and maximum, and values normalized by the growth rate. Reactions are sorted by FBA flux value.

#### 7. Fluxes per pathway

For each pathway specified in `pathwaysOfInterest`, the report includes:

- **Fluxes** — FBA, FVA min/max, and normalized values for all reactions in the pathway
- **FBA flux sum** — aggregate flux sum computed from the FBA solution
- **Sampling flux sum** — aggregate flux sum computed from sampling results (if sampling was performed)

!!! note "Pathway matching"
    Pathway names must match the `subSystems` field of the model. Pathways not found in the model are skipped with a console warning.

#### 8. Fluxes per metabolite

For each metabolite specified in `metsOfInterest`, the report includes:

- **Fluxes** — FBA, FVA min/max, and normalized values for all reactions involving the metabolite
- **FBA flux sum** — aggregate flux sum computed from the FBA solution
- **Sampling flux sum** — aggregate flux sum computed from sampling results (if sampling was performed)

!!! note "Metabolite matching"
    Metabolite IDs are matched using a prefix pattern (e.g. `glc_D` matches `glc_D[e]`). Metabolites not found in the model are skipped with a console warning.

### Output file

The PDF file is named `<modelName>_AnalysisReport_<analysisId>.pdf` and saved to the directory specified by `path`. If the directory does not exist, it is created automatically.

### Usage example

```matlab
% Generate a report with default settings
writeAnalysisReport(project, "model1", "analysis_20240815_1430", ...
    'path', './results/reports/');

% Generate a report focusing on specific pathways and metabolites
writeAnalysisReport(project, "model1", "analysis_20240815_1430", ...
    'path', './results/reports/', ...
    'pathwaysOfInterest', {"Glycolysis", "TCA cycle"}, ...
    'metsOfInterest', {"glc_D", "o2", "ac"});
```

!!! note "Upcoming features"
    The following report sections are planned for future integration:

    - **Shadow prices** — metabolite shadow price analysis
    - **Figures and plots** — embedded visualizations of flux distributions and sampling landscapes