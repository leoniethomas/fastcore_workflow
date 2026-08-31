# Input data format

To create a `project`, a `params` structure must be provided as input to the pipeline. The table below describes the expected format for each field.

## Required fields

| Field | Type | Format | Description |
|-------|------|--------|-------------|
| `modelName` | char / string | scalar | Name of the model. Becomes the key in `project.models`. |
| `contextSpecificModel` | struct (COBRA) | scalar | Genome-scale metabolic model in COBRA format (`.S`, `.rxns`, `.genes`, `.rules`, `.lb`, `.ub`, etc.). |

## Optional fields

| Field | Type | Format | Description |
|-------|------|--------|-------------|
| `objFunction` | char / string | scalar | Name of the objective reaction (e.g. biomass). |
| `referenceModel` | char / string | scalar | Name of the reference model for comparisons. |
| `consensusProportion` | numeric | scalar, ≥ 0 | Proportion of samples that must express a gene for it to be included in a consensus model (default: 0.9). |
| `mediumComposition` | table | — | Medium composition (metabolites and exchange bounds). |
| `mapping` | double | matrix (sparse or dense) | Gene-to-reaction mapping matrix. |
| `coreReactions` | cell array | vector or empty | Cell array of core reaction identifiers. |
| `optionalSettings` | struct | scalar | Additional settings. See [optional settings](#optional-settings) below. |

## Conditional fields

These fields are required only when expression or discretized data is provided.

| Field | Type | Format | Required when | Description |
|-------|------|--------|---------------|-------------|
| `discretizedData` | double | 2-D matrix (genes × samples) | Optional, but if present triggers dependencies | Discretized expression values (−1, 0, 1). |
| `expressionData` | double | 2-D matrix (genes × samples) | Optional, but if present triggers dependencies | Raw expression values (e.g. TPM or FPKM). |
| `geneIds` | string array | column vector (char / cellstr accepted and auto-converted) | **Required** if `discretizedData` or `expressionData` is present | Gene identifiers matching the rows of the expression data. |
| `dico` | table | must contain columns `geneIdsInModel` and `geneIdsInData` | **Required** if `discretizedData` or `expressionData` is present | Mapping table between data gene IDs and model gene IDs. |
| `sampleMetadata` | table | — | Must be provided **together with** `sampleLabeling` (both or neither) | Sample metadata (conditions, groups, etc.). |
| `sampleLabeling` | char / string | scalar, must match a column name in `sampleMetadata` | Must be provided **together with** `sampleMetadata` (both or neither) | Name of the column in `sampleMetadata` that defines sample groups. |
| `manuallySetBoundaries` | struct | scalar | Optional | Manually defined exchange boundaries. See [manually set boundaries](#manually-set-boundaries) below. |

## Cross-field validation rules

| Rule | Condition |
|------|-----------|
| `sampleMetadata` ↔ `sampleLabeling` | Must be provided together or both omitted. `sampleLabeling` must match a column name in `sampleMetadata`. |
| `geneIds` ↔ `discretizedData` | `geneIds` must have the same number of rows as `discretizedData`. |
| `geneIds` ↔ `expressionData` | `geneIds` must have the same number of rows as `expressionData`. |
| `discretizedData` ↔ `expressionData` | If both present, must have the same dimensions (rows = genes, columns = samples). |
| `dico.geneIdsInModel` ↔ `contextSpecificModel.genes` | `dico.geneIdsInModel` must cover ≥ 5% of `contextSpecificModel.genes`. A warning is issued with the actual coverage percentage. |
| `dico.geneIdsInData` ↔ `geneIds` | `dico.geneIdsInData` must cover ≥ 5% of `geneIds`. A warning is issued with the actual coverage percentage. |

## Sub-structures

<a id="optional-settings"></a>

### `optionalSettings` (struct)

| Field | Type | Description |
|-------|------|-------------|
| `medium` | — | Optional medium definition. |
| `notMediumConstrained` | — | Reactions not constrained by the medium. |
| `func` | — | Reaction(s) forced to carry a flux. |

<a id="manually-set-boundaries"></a>

### `manuallySetBoundaries` (struct)

| Field | Type | Description |
|-------|------|-------------|
| `closedImports` | — | Import reactions to close. |
| `closedExports` | — | Export reactions to close. |
| `unconstrainedImports` | — | Import reactions to leave unconstrained. |
| `unconstrainedExports` | — | Export reactions to leave unconstrained. |

### `dico` (table)

| Column | Type | Description |
|--------|------|-------------|
| `geneIdsInModel` | string / char / cellstr | Gene identifiers as they appear in `contextSpecificModel.genes`. Must cover ≥ 5% of model genes. |
| `geneIdsInData` | string / char / cellstr | Gene identifiers as they appear in `geneIds` / expression data. Must cover ≥ 5% of data gene IDs. |

!!! tip "Building the dico table"
    The `dico` table can be assembled using [Ensembl BioMart](https://www.ensembl.org/biomart/martview).
    It maps the gene identifiers from your expression data to the gene identifiers
    in the COBRA model, which may use different naming conventions.

## Minimal example

The minimal input requires only a model name and a COBRA model:

```matlab
paramsForPipeline.modelName = 'myModel';
paramsForPipeline.contextSpecificModel = myCobraModel;
```

A typical input with expression data:

```matlab
paramsForPipeline.modelName = 'myModel';
paramsForPipeline.contextSpecificModel = myCobraModel;
paramsForPipeline.expressionData = tpmMatrix;       % genes × samples (double)
paramsForPipeline.geneIds = geneIdStrings;          % string array, one per row
paramsForPipeline.dico = dicoTable;                 % table with geneIdsInModel, geneIdsInData
paramsForPipeline.sampleMetadata = sampleTable;     % sample info
paramsForPipeline.sampleLabeling = 'condition';     % column name in sampleMetadata
paramsForPipeline.objFunction = 'biomass_reaction';
paramsForPipeline.consensusProportion = 0.9;

project = createProject(paramsForPipeline);
```