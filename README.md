# ShinyCellModular
***
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

**Note:** [FlexDotPlot](https://github.com/Simon-Leonard/FlexDotPlot) is required but not on CRAN. Install it with:
```r
devtools::install_github("Simon-Leonard/FlexDotPlot")
```

**Note:** `limma` and `edgeR` are Bioconductor packages. Install them with:
```r
 BiocManager::install(c("limma", "edgeR"))
```

**Note:** [ShinyCell](https://github.com/SGDDNB/ShinyCell) is required but not on CRAN. Install it with:
```r
devtools::install_github("SGDDNB/ShinyCell")
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
                     )
```

### 3. `useShinyCellModular()`

```r
# Create a new app.R with the modular ShinyCellModular tabs

useShinyCellModular(
    shiny.dir = "testing_data/",
    rsconnect.deploy = FALSE,
    data_type = "RNA",
    overwrite_modules = TRUE,
    app_title = "Testing"
)

runApp("testing_data")
# or open app.R and run
```

***
## Legacy version

The pre-package version of ShinyCellModular is preserved in the [`legacy` branch](https://github.com/MonashBioinformaticsPlatform/ShinyCellModular/tree/legacy) for users who are already working with that code. New development happens on `main`.

***
## Acknowledgement

We would love to know if ShinyCellModular is useful to you and your team. If you use it in your work or build new modules on top of it, please let us know and acknowledge it in your publications — this helps us track its impact and justify continued development.


