# Installation

This page describes the software requirements and installation steps for the GEM-PATCH analysis pipeline.

## MATLAB version

The pipeline requires **MATLAB R2024a or R2024b**.

!!! warning "MATLAB R2025 is not compatible"
    The UMAP toolbox used by this pipeline relies on Java access to MATLAB figures, which is removed in MATLAB R2025. Use R2024b or earlier.

Check your MATLAB version in the command window:

```matlab
ver
```

## Required toolboxes

The following MathWorks toolboxes must be installed. Install them via **Home > Add-Ons > Get Add-Ons** in MATLAB, or through your institution's license portal.

| Toolbox | Required for | Minimum version |
|---|---|---|
| Statistics and Machine Learning Toolbox | UMAP, statistical analysis | R2024a |
| MATLAB Report Generator | Automated report generation | R2024a |

Verify that all toolboxes are installed:

```matlab
ver
```

The output should list both `Statistics and Machine Learning Toolbox` and `MATLAB Report Generator`.

## UMAP toolbox installation

The pipeline uses the MATLAB UMAP implementation by Meehan et al. [[File Exchange #71902](https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)].

### Step 1 — Download

Download the toolbox from [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap).

Alternatively, open it directly in MATLAB Online by clicking **Open in MATLAB Online** on the File Exchange page.

### Step 2 — Add to MATLAB path

Unzip the downloaded archive to a permanent location (e.g. `C:\Tools\umap`), then add it to your MATLAB path:

```matlab
addpath(genpath('C:\Tools\umap'));
savepath;
```

!!! tip "Save the path"
    Run `savepath` after `addpath` so the UMAP toolbox is available in future MATLAB sessions without re-adding it manually.


### Step 3 — Verify installation

```matlab
% Should run without errors and produce a 2D scatter plot
X = rand(100, 10);
Y = run_umap(X, 'n_components', 2);
scatter(Y(:,1), Y(:,2));
```

!!! note "Dependencies"
    The UMAP toolbox requires the **Statistics and Machine Learning Toolbox** and the **Financial Toolbox**. Check with `ver` if you encounter missing function errors.

## Report generation

The pipeline generates automated analysis reports using the **MATLAB Report Generator** toolbox. Reports are produced by the `writeAnalysisReport` function.

### How it works

The report generator uses the Report API and Document Object Model (DOM) API to programmatically build documents from your analysis results. The pipeline captures tables, and text output from the single model analysis step and compiles them into a structured PDF report. The report file is saved to the working directory.

## Verification checklist

After completing all installation steps, verify your setup with the following commands in the MATLAB command window:

```matlab
% 1. Check MATLAB version (should be R2024a or R2024b)
ver('matlab')

% 2. Check required toolboxes
ver('stats')      % Statistics and Machine Learning Toolbox
ver('rptgen')     % MATLAB Report Generator

% 3. Check UMAP installation
which run_umap    % Should return the path to run_umap.m

% 4. Test UMAP
Y = run_umap(rand(50, 5), 'n_components', 2);
disp(size(Y))     % Should be [50, 2]
```

If all commands execute without errors, your environment is ready for the pipeline.

## Building models with rFastcormics

rFASTCORMICS will be integrated into the next release of the COBRA Toolbox. In the meantime, if you want to build context-specific models using rFASTCORMICS before comparing them with GEM-PATCH, follow the installation instructions on the [rFASTCORMICS GitHub page](https://github.com/sysbiolux/rFASTCORMICS/tree/master/rFASTCORMICS%20for%20RNA-seq%20data/rFASTCORMICS_v2).

A supplementary `rFastcormicsPipeline.m` script is provided in the `scr/` folder. It builds context-specific models and outputs a `params` structure that can be passed directly to `createProject.m`. For a usage example, see `no1_modelBuildingAndProjectInit.m` in the `BRCAexample` folder.