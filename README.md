# ShinyCellPlus
***
**ShinyCellPlus** is a modular version of [ShinyCell](https://github.com/SGDDNB/ShinyCell) developed at the Monash Genomics and Bioinformatics Platform (MGBP). Each module consists on a tab in the app. Each module is created individually and is it selfcontained. **ShinyCellPlus** supports large scRNAseq and multimodal datasets with fast on-demand HDF5 access, extended visualisations, improved filtering, and publication-ready plots. Its modular structure makes it flexible, scalable, and easy to customise.

## Features

- Modular UI and server structure
- Supports scRNAseq, ATAC, and multimodal datasets
- Fast HDF5 on-demand loading
- Publication‑ready plots (PNG/PDF export)
- Extended visualisation tabs (UMAP, 3D UMAP, violin, bubble, heatmap, coexpression, AUC marker genes)
- Cell subsetting and conditional plotting
- Easy integration with new modules via a registry system

***
## Fast usage just needs 3 steps

### 1. Setup

Clone this repository

```
git clone https://github.com/MonashBioinformaticsPlatform/ShinyCellPlus.git
```

Open the **.Rproj** file

Load RENV - all require library

```
install.packages("renv")
renv::restore()
```


Run the 2 helper functions `prepShinyCellPlus()` and `useShinyCellPlus()`



### 2. `prepShinyCellPlus()`

```
library(ShinyCell) #devtools::install_github("SGDDNB/ShinyCell")
library(Seurat)
library(Signac)
library(dplyr)


# Prepare seurat object, checks Key names, creates sc1counts.h5, adds a 3D reduction UMAP
cnts<-LoadSeuratRds("seurat_object.Rds")

source("functions/prepShinyCellPlus.R")

prepShinyCellPlus(seurat_rds = "seurat_obj.rds", # or seurat_obj = cnts
                  out_dir = "testing_data", 
                  do_umap3d = TRUE,  
                  do_markers= TRUE,   
                  markers_res_pattern = "RNA_snn_res")

```


### 3.`useShinyCellPlus()`

```

# Create a new app.R with the modified ShinyCell tabs

source("functions/useShinyCellPlus.R")

useShinyCellPlus(
    shiny.dir="testing_data/",
    shinycellplus.dir.src="~/ShinyCellPlus/",
    rsconnect.deploy = FALSE,
    data_type = "",
    enabled_tabs = c("cellinfo_cellinfo",
                    "cellinfo_geneexpr",
                    "cellinfo3D_cellinfo3D",
                    "cellinfo3D_geneexpr3D",
                    "genecoex",
                    "violin_boxplot",
                    "proportions",
                    "bubble_heatmap",
                    "pseudobulk"),
    overwrite_modules = TRUE,
    app_title='Testing'
)


```
Review Docs for further information on development instructions 