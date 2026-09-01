# Single Model Analysis

The `singleModelAnalysis` function runs one or more analyses on a model or a list of models and stores the results in the project structure.

## Prerequisites

!!! warning "Project format"
    The `project` structure must follow the format described in the [Project Layout](project_layout.md) page. Use `createProject` and `addModelsToProject` to build a valid project before running analyses.

## Available analyses

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

## Function signature

```matlab
project = singleModelAnalysis(project, parameterTable, modelList, analyses, saveCheckpoint, resumeFromCheckpoint)
```

### Input arguments

| Name | Type | Description | Default |
|---|---|---|---|
| `project` | `struct` | Project structure following the format described in [Project Layout](project_layout.md) | required |
| `parameterTable` | `table` | Parameter table defining analysis settings (at minimum the default one) | required |
| `modelList` | `string array` | Names of the models to analyze | `{}` (all models) |
| `analyses` | `cell array` | List of analyses to perform | `{}` (all available) |
| `saveCheckpoint` | `logical` | Whether to save a checkpoint after each model | `true` |
| `resumeFromCheckpoint` | `logical` | Whether to resume from the last saved checkpoint | `false` |

### Output

| Name | Type | Description |
|---|---|---|
| `project` | `struct` | The input project with an `analysis` field added to each model |

## Analysis IDs

Each analysis run is assigned a unique identifier based on the current timestamp:

```
analysis_20240815_1430
```

The format is `analysis_yyyyMMdd_HHmm`, where:

- `yyyyMMdd` — date (year, month, day)
- `HHmm` — time (hour, minute)

This ID is used as a field name under `project.models.<modelName>.analysis.<id>` and allows multiple analysis runs to coexist on the same model without overwriting previous results. When comparing models later, you can reference a specific analysis by its ID.

## Parameter storage

For each analysis run, the parameter table used is stored alongside the results:

```matlab
project.models.<modelName>.analysis.<id>.parameters
```

This ensures full traceability — you can always retrieve the exact settings that produced a given set of results, even if the parameter table was modified between runs.

## Loopless sampling

Loopless sampling can be run in two ways:

1. **In the same run as sampling** — if both `sampling` and `loopless` are requested, the sampling results from the current run are automatically used as input for the loopless analysis.

2. **Using a previous sampling run** — if only `loopless` is requested (without `sampling`), you must specify which existing sampling results to use via the `samplingToUse` parameter in the parameter table. This parameter should contain the analysis ID of a previously run sampling analysis, in the same order as the models in `modelList`.

!!! warning "Order matters"
    When using `samplingToUse`, the analysis IDs must be listed in the same order as their corresponding models in `modelList`. Mismatched ordering will produce an error.

## Checkpoint system

The function includes a checkpoint mechanism to protect against crashes during long-running analyses:

- **Automatic save** — after each model completes, the current project state and progress index are saved to `singleModelAnalysis_checkpoint.mat` (if `saveCheckpoint` is `true`).
- **Resume** — set `resumeFromCheckpoint = true` to load the last checkpoint and continue from the next model. The function prints which model it is resuming from.
- **Error handling** — if an error occurs during analysis of a model, the function catches it, prints the error message and stack trace, and returns the partial project. You can then relaunch with `resumeFromCheckpoint = true` to continue.
- **Cleanup** — the checkpoint file is automatically deleted once all models have been successfully analyzed.

## Usage example

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