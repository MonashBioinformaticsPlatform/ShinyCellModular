<p align="center">
  <img src="images/ShinyCellModularLogo_banner.png" width="800">
</p>

***
## Motivation

Single-cell apps rarely stay static once a lab starts using them: labels get reordered, axes get fixed, colours get matched to a figure, a new tab gets requested. If you're running a ShinyCell-based app for more than one project, this will sound familiar: each of those small changes usually means touching a large, interconnected script. That is what **ShinyCellModular** is for: an architecture that lets us keep tailoring single-cell apps to individual researchers and projects, quickly and without waiting on upstream changes.

Most ShinyCell-derived apps are one large, cross-referenced script: tabs share state and call into each other, so touching one thing risks breaking another. ShinyCellModular restructures this as a plugin system: every tab is a self-contained R file (UI + server + `register_tab()`), discovered automatically at build time from a directory listing, with no hard-coded catalogue to edit. Add a tab by dropping in a file; remove one without touching anything else; hand a tab to a team member without them needing to understand the rest of the app.

That architecture is what lets us keep customising ShinyCellModular for our own researchers over time, rather than accumulating one-off patches that get harder to maintain with every request. We build most tabs ourselves, driven by real project needs, but the module structure is deliberately open: see [creating your own modules/tabs](documentation/developer_guide.md) if you want to build on it or adapt it for your own group.

## What it does? 

**ShinyCellModular** is an R package, a modular version of [ShinyCell](https://github.com/SGDDNB/ShinyCell) developed at the Monash Genomics and Bioinformatics Platform (MGBP). It takes your [Seurat](https://satijalab.org/seurat/) object from single cell experiments and creates an interactive Shiny app to explore your data. Each module is a tab in the app, created individually and self-contained. **ShinyCellModular** supports large scRNAseq and multimodal datasets with fast on-demand HDF5 and parquet access, extended visualisations, improved filtering, and publication-ready plots. Its modular structure makes it flexible, scalable, and easy to customise and to patch.

[Example of ShinyCellModular app and tutorials](https://bioinformatics.erc.monash.edu/shinyapps-public/app/scrnaseq-shinycellmodular-example)
 
Review Docs for further information on [developer guide](documentation/developer_guide.md)     

## Features

- Modular UI and server structure
- Supports scRNAseq, ATAC, and multimodal datasets
- Fast HDF5 and parquet on-demand loading
- Publication-ready plots (PNG/PDF export)
- Extended visualisation tabs (UMAP, 3D UMAP, violin, bubble, heatmap, coexpression, marker genes)
- Pseudobulk differential expression
- Cell subsetting and conditional plotting
- Marker gene visualisation from precomputed parquet files
- Per-tab authorship and metadata footer
- Easy integration with new modules via a registry system
- Deployment to Posit Connect via rsconnect

***
## Fast usage just needs 3 steps

### 1. Setup

Install the package directly from GitHub:

```r
devtools::install_github("MonashBioinformaticsPlatform/ShinyCellModular")
library(ShinyCellModular)
```

The first time you run `prepShinyCellModular` add `install_missing = TRUE` to auto-install any missing dependencies:

```r
prepShinyCellModular(install_missing = TRUE)
```

Run the 2 helper functions `prepShinyCellModular()` and `useShinyCellModular()`

### 2. [`prepShinyCellModular()`](documentation/prepShinyCellModular_explained.md)

```r
library(ShinyCellModular)

# Prepare seurat object, checks Key names, creates sc1counts.h5, adds a 3D UMAP reduction, identify marker genes for all resolutions
prepShinyCellModular(seurat_rds = "seurat_object.rds", # or seurat_obj = cnts,
                     out_dir = "testing_data_RNA", 
                     assays_selected = "RNA", 
                     do_umap3d = TRUE,  
                     do_markers = TRUE
                     #, install_missing = TRUE
                     )
```

### 3. [`useShinyCellModular()`](documentation/useShinyCellModular_explained.md)

```r
# Create a new app.R with the modular ShinyCellModular tabs

useShinyCellModular(
    out_dir  = "testing_data/",
    data_type = "RNA",
    overwrite_modules = TRUE, # be careful with this if you have done any manual changes to the modules code, it will replace the whole folder with the package modules code
    app_title = "Testing"
)

runApp("testing_data")
# or open app.R and run
```

To include only specific tabs pass their IDs to `enabled_tabs`:

```r
useShinyCellModular(
                    shiny.dir    = "testing_data/",
                    data_type    = "RNA",
                    enabled_tabs = c("cellinfo_cellinfo", "violin_boxplot", "pseudobulk"),
                    app_title    = "Testing"
)
```

***
## Available tabs

### RNA

| Tab ID (`enabled_tabs`)   | Tab title              | What it shows                                         | Extra prep needed             |
|---------------------------|------------------------|-------------------------------------------------------|-------------------------------|
| `cellinfo_cellinfo`       | CellInfo vs CellInfo   | 2D embedding coloured by metadata                     | n/a                          |
| `cellinfo_geneexpr`       | CellInfo vs GeneExpr   | 2D embedding with gene expression overlay             | n/a                          |
| `cellinfo3D_cellinfo3D`   | CellInfo3D             | Interactive 3D embedding coloured by metadata         | `do_umap3d = TRUE` in prep    |
| `cellinfo3D_geneexpr3D`   | CellInfo3D vs GeneExpr | Interactive 3D embedding with gene expression overlay | `do_umap3d = TRUE` in prep    |
| `genecoex`                | Gene Coexpression      | Coexpression of selected genes across cells or groups | n/a                          |
| `violin_boxplot`          | Violin / BoxPlot       | Violin and boxplots for gene expression or metadata   | n/a                          |
| `proportions`             | Cell Proportions       | Cell composition across groups                        | n/a                          |
| `bubble_heatmap`          | Bubble Plot / Heatmap  | Bubble plot and heatmap for gene sets across groups   | n/a                          |
| `pseudobulk`              | Pseudobulk DE          | Pseudobulk aggregation and differential expression    | `do_counts_h5 = TRUE` in prep |

### QC function
Coming soon.

### ATAC

Coming soon.

### RNA_ATAC

Coming soon.

### SPATIAL

Coming soon.

### CropSeq/ PerturbSeq

Coming soon.


***
## Legacy version

The pre-package version of ShinyCellModular is preserved in the [`legacy` branch](https://github.com/MonashBioinformaticsPlatform/ShinyCellModular/tree/legacy) for users who are already working with that code. New development happens on `main`.

***
## Acknowledgement

We'd love to hear if ShinyCellModular is useful to you outside MGBP. If you use it in your work or build new modules on top of it, please let us know and acknowledge it in your publications; this helps us track its impact and justify continued development.
