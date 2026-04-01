# ShinyCellModular
***
**ShinyCellModular** is a modular version of [ShinyCell](https://github.com/SGDDNB/ShinyCell) developed at the Monash Genomics and Bioinformatics Platform (MGBP). Each module consists on a tab in the app. Each module is created individually and is it selfcontained. **ShinyCellModular** supports large scRNAseq and multimodal datasets with fast on-demand HDF5 access, extended visualisations, improved filtering, and publication-ready plots. Its modular structure makes it flexible, scalable, and easy to customise. [Example of ShinyCellModular in Single Cell RNAseq](https://bioinformatics3.erc.monash.edu/rsconnect/content/543/)

Review Docs for further information on [functions details](docs/functions_details.md)    
Review Docs for further information on [development instructions](docs/developer_guide.md)     

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
git clone https://github.com/MonashBioinformaticsPlatform/ShinyCellModular.git
```

Open the **.Rproj** file

Load RENV - all require library

```
install.packages("renv")
renv::restore()
```


Run the 2 helper functions `prepShinyCellModular()` and `useShinyCellModular()`



### 2. `prepShinyCellModular()`

```
library(ShinyCell) #devtools::install_github("SGDDNB/ShinyCell")
library(Seurat)
library(Signac)
library(dplyr)


# Prepare seurat object, checks Key names, creates sc1counts.h5, adds a 3D reduction UMAP
cnts<-LoadSeuratRds("seurat_object.Rds")

source("functions/prepShinyCellModular.R")

prepShinyCellModular(seurat_rds = "seurat_object.rds",
                  out_dir = "testing_data_RNA", 
                  assays_selected = "RNA",
                  do_umap3d = TRUE,  
                  do_markers= TRUE,   
                  markers_res_pattern = "wsnn_res",
                  )


```


### 3.`useShinyCellModular()`

```

# Create a new app.R with the modified ShinyCell tabs

source("functions/useShinyCellModular.R")

useShinyCellModular(
    shiny.dir="testing_data/",
    shinycellmodular.dir.src="~/ShinyCellModular/",
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


