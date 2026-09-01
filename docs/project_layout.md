# Project layout

A `project` is a MATLAB hierarchical structure. Click on each field to discover the associated subfields.

!!! warning "Field format differences"
    Field formats may differ between what is required to initialize a project and how they are stored downstream in the `project`. Before building the pipeline, check the [Project Initialization](project_init.md) page for field requirements.

???+ note "project (struct)"

    | Field | Type | Status | Description |
    |-------|------|--------|-------------|
    | `models` | struct | **Required** | One or more models |
    | `comparisons` | struct | Optional | Comparisons between models |

??? note "project.models (struct)"

    Each field is a model name (`<modelName>`).

    ???+ note "project.models.&lt;modelName&gt; (struct)"

        | Field | Type | Status | Dimensions | Description |
        |-------|------|--------|------------|-------------|
        | `model` | struct (COBRA) | **Required** | — | Metabolic model |
        | `sampleMetadata` | table | Optional | nbSamples × N | Sample metadata |
        | `discretizedData` | table | Optional | nbGenes × 2+ | Discretized data (`geneIds`, `value`) |
        | `expressionData` | table | Optional | nbGenes × 2+ | Raw expression data (`geneIds`, `expression`) |
        | `mappedDiscretizedRxnsAllSamples` | int8 | Optional | nbRxns × N | Mapping across all samples |
        | `mappedDiscretizedRxns` | int8 | Optional | nbRxns × N | Consensus mapping |
        | `coreReactions` | vector | Optional | — | Core reaction indices |
        | `settings` | struct | Optional | — | Configuration parameters |
        | `analysis` | struct | Optional | — | Analysis results |

        ??? note "settings (struct)"

            | Field | Type | Status | Description |
            |-------|------|--------|-------------|
            | `medium` | struct | Optional | Medium composition |
            | `scriptParameters` | struct | Optional | Script parameters |
            | `dico` | table | Conditional* | Gene ID mapping (data ↔ model) |
            | `objFunction` | char | Optional | Objective reaction |
            | `referenceModel` | char | Optional | Reference model |
            | `mapping` | sparse double | Optional | Mapping matrix |
            | `optionalSettings` | struct | Optional | Additional parameters |

            *\*Required if `discretizedData` or `expressionData` is present.*

            ??? note "settings.dico (table)"

                | Column | Status | Description |
                |---------|--------|-------------|
                | `geneIdsInModel` | **Required** | Model gene IDs (same order as `model.genes`) |
                | `geneIdsInData` | Conditional | Data gene IDs (same order as `discretizedData.geneIds`) |

            ??? note "settings.medium (struct)"

                | Field | Type | Description |
                |-------|------|-------------|
                | `mediumComposition` | table | Medium composition |
                | `manuallySetBoundaries` | struct | Manually set boundaries |

                ??? note "medium.manuallySetBoundaries (struct)"

                    | Field | Description |
                    |-------|-------------|
                    | `closedImports` | Closed imports |
                    | `closedExports` | Closed exports |
                    | `unconstrainedImports` | Unconstrained imports |
                    | `unconstrainedExports` | Unconstrained exports |

            ??? note "settings.scriptParameters (struct)"

                | Field | Type | Description |
                |-------|------|-------------|
                | `consensusProportion` | numeric | Consensus proportion (default: 0.9) |
                | `sampleLabeling` | string | Column in `sampleMetadata` defining sample groups |

            ??? note "settings.optionalSettings (struct)"

                | Field | Description |
                |-------|-------------|
                | `medium` | Optional medium |
                | `notMediumConstrained` | Reactions not constrained by medium |
                | `func` | Forced reaction(s) |

        ??? note "analysis (struct)"

            One entry per analysis (`analysis_&lt;id&gt;`) plus an `active` entry.

            ???+ note "analysis.&lt;id&gt; (struct)"

                | Field | Type | Status | Description |
                |-------|------|--------|-------------|
                | `parameters` | table | Optional | Parameters (`Parameter`, `Analysis`, `Value`) |
                | `FBA` | struct | Optional | Flux Balance Analysis |
                | `FVA` | struct | Optional | Flux Variability Analysis |
                | `singleGeneDeletion` | struct | Optional | Single gene deletion |
                | `doubleGeneDeletion` | struct | Optional | Double gene deletion |
                | `sampling` | struct | Optional | Flux sampling |
                | `kld` | struct | Optional | Kullback-Leibler Divergence |

                !!! tip "`analysisId` field"
                    The `active` entry contains an additional `analysisId` (char) 
                    field in each sub-struct.

                ??? note "FVA (struct)"

                    | Field | Type | Dimensions | Description |
                    |-------|------|------------|-------------|
                    | `minMaxFluxes` | table | nbRxns × 2 | Columns `minFlux`, `maxFlux` |
                    | `loopStatus` | logical | nbRxns × 1 | Loop status |
                    | `analysisId` | char | — | *(active entry only)* |

                ??? note "singleGeneDeletion (struct)"

                    | Field | Description |
                    |-------|-------------|
                    | `grRatio` | Growth ratio KO/WT |
                    | `grRateKO` | Growth rate after deletion |
                    | `grRateWT` | Wild-type growth rate |
                    | `hasEffect` | Whether deletion has an effect |
                    | `delRxns` | Deleted reactions |
                    | `fluxSolution` | Flux solution |

                ??? note "doubleGeneDeletion (struct)"

                    | Field | Description |
                    |-------|-------------|
                    | `grRatioDble` | Double deletion growth ratio |
                    | `grRatioKO` | Growth ratio after double deletion |
                    | `grRateWT` | Wild-type growth rate |

                ??? note "sampling (struct)"

                    | Field | Type | Dimensions | Description |
                    |-------|------|------------|-------------|
                    | `modelSampling` | struct | — | Sampling model |
                    | `samples` | single | nbRxns × nSamp | Flux samples |
                    | `cycleFreeFlux` | struct | — | Cycle-free results |

                    ??? note "sampling.cycleFreeFlux (struct)"

                        | Field | Type | Dimensions | Description |
                        |-------|------|------------|-------------|
                        | `samplesLl` | single | nbRxns × nSamp | Loopless samples |
                        | `thermoFeas` | uint8 | nbRxns × nSamp | Thermodynamic feasibility |
                        | `sampleStatusAfterCorrection` | uint8 | nbRxns × nSamp | Status after correction |
                        | `neededAttempts` | uint8 | nbRxns × 1 | Number of attempts needed |
                        | `looplessStatus` | uint8 | nbRxns × nSamp | Loopless status |

                ??? note "kld (struct)"

                    | Field | Description |
                    |-------|-------------|
                    | `samplingSets` | Compared sample sets |
                    | `kldMatrix` | KLD matrix |
                    | `pValueKld` | p-value |
                    | `fdr` | False Discovery Rate |
                    | `setLabels` | Set labels |

??? note "project.comparisons (struct)"

    Each field is named `Name1_vs_Name2__date`.

    ???+ note "comparisons.&lt;Name1_vs_Name2__date&gt; (struct)"

        | Field | Type | Dimensions | Description |
        |-------|------|------------|-------------|
        | `modelNames` | string | nbModels × 1 | Compared models |
        | `referenceModel` | string | — | Reference model |
        | `analysisIds` | table | nbAnalyses × nbModels | Compared analyses |
        | `structuralComparison` | struct | — | Structural comparison |
        | `functionalComparison` | struct | — | Functional comparison |
        | `samplingComparison` | struct | — | Sampling comparison |

        ??? note "structuralComparison (struct)"

            | Field | Type | Dimensions | Description |
            |-------|------|------------|-------------|
            | `rxnMappingTable` | table | nbRxns × nbModels | Reaction mapping |
            | `plots` | struct | — | Generated plots |

            ??? note "structuralComparison.plots"

                | Field | Description |
                |-------|-------------|
                | `dataDiscretization` | Data discretization |
                | `coreReactions` | Core reactions |
                | `coreReactionsIntersections` | Core reaction intersections |
                | `intersections` | Global intersections |
                | `jaccardDist` | Jaccard distance |
                | `reactionPathwayPresence` | Reaction presence by pathway |

        ??? note "functionalComparison (struct)"

            | Field | Description |
            |-------|-------------|
            | `plots.objValue` | Objective value |
            | `plots.import` | Import fluxes |
            | `plots.export` | Export fluxes |
            | `plots.fvaSim.overall` | Overall FVA similarity |
            | `plots.fvaSim.hist` | FVA similarity histogram |
            | `plots.fvaSim.enrich` | FVA similarity enrichment |
            | `plots.fba.heatmapRxnFluxsum` | Heatmap of flux sum per reaction |
            | `plots.fba.heatmapRxnActivity` | Heatmap of reaction activity |
            | `plots.fba.heatmapMetsFluxsum` | Heatmap of flux sum per metabolite |

        ??? note "samplingComparison (struct)"

            | Field | Type | Dimensions | Description |
            |-------|------|------------|-------------|
            | `orderedSamples` | double | nbRxns × cumSamp | Ordered samples |
            | `sampleModelLabels` | string | 1 × cumSamp | Model label per sample |
            | `orderedFba` | double | nbRxns × nbModels | Ordered FBA solutions |
            | `plots` | struct | — | Generated plots |

            ??? note "samplingComparison.plots"

                | Field | Description |
                |-------|-------------|
                | `heatmapRxnFluxSum` | Heatmap of flux sum per reaction |
                | `heatmapMetsFluxSum` | Heatmap of flux sum per metabolite |
                | `heatmapRxnFluxSumSamples` | Heatmap of flux sum per reaction (samples) |
                | `heatmapMetsFluxSumSamples` | Heatmap of flux sum per metabolite (samples) |