# ShinyCellModular
***

**ShinyCellModular** takes your [Seurat](https://satijalab.org/seurat/) object from single cell experiments and creates an interactive Shiny app to explore your data.

---

**ShinyCellModular** is an R package, a modular version of [ShinyCell](https://github.com/SGDDNB/ShinyCell) developed at the Monash Genomics and Bioinformatics Platform (MGBP). Each module is a tab in the app, created individually and self-contained. **ShinyCellModular** supports large scRNAseq and multimodal datasets with fast on-demand HDF5 and parquet access, extended visualisations, improved filtering, and publication-ready plots. Its modular structure makes it flexible, scalable, and easy to customise and to patch.

[Example of ShinyCellModular app and tutorials](https://bioinformatics3.erc.monash.edu/rsconnect/content/543/)

Review Docs for further information on [functions details](docs/functions_details.md)    
Review Docs for further information on [development instructions](docs/developer_guide.md)     

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

### 2. `prepShinyCellModular()`

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

### 3. `useShinyCellModular()`

```r
# Create a new app.R with the modular ShinyCellModular tabs

useShinyCellModular(
    shiny.dir = "testing_data/",
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
| `cellinfo_cellinfo`       | CellInfo vs CellInfo   | 2D embedding coloured by metadata                     | —                             |
| `cellinfo_geneexpr`       | CellInfo vs GeneExpr   | 2D embedding with gene expression overlay             | —                             |
| `cellinfo3D_cellinfo3D`   | CellInfo3D             | Interactive 3D embedding coloured by metadata         | `do_umap3d = TRUE` in prep    |
| `cellinfo3D_geneexpr3D`   | CellInfo3D vs GeneExpr | Interactive 3D embedding with gene expression overlay | `do_umap3d = TRUE` in prep    |
| `genecoex`                | Gene Coexpression      | Coexpression of selected genes across cells or groups | —                             |
| `violin_boxplot`          | Violin / BoxPlot       | Violin and boxplots for gene expression or metadata   | —                             |
| `proportions`             | Cell Proportions       | Cell composition across groups                        | —                             |
| `bubble_heatmap`          | Bubble Plot / Heatmap  | Bubble plot and heatmap for gene sets across groups   | —                             |
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

We would love to know if ShinyCellModular is useful to you and your team. If you use it in your work or build new modules on top of it, please let us know and acknowledge it in your publications — this helps us track its impact and justify continued development.
