
# Developer Instructions


## prepShinyCellPlus()
  
  `prepShinyCellPlus()` is the preparation step that builds everything the ShinyCellPlus app will need later. It takes a Seurat object, optionally adds extra artefacts (markers, 3D UMAP, counts H5), and then runs ShinyCell to generate the standard ShinyCell output folder. The result is a ready to use directory on disk (default `Files_ShinyCell`) containing the `.rds` files and optional extras that `useShinyCellPlus()` and the modular tabs can consume.

Step by step, what it does:
  
  1. Sets defaults for output files
If you do not provide paths:
  
  * `markers_file` defaults to `out_dir/markergenes_lists.parquet`
* `counts_h5_file` defaults to `out_dir/sc1counts.h5`

2. Checks dependencies and optionally installs them
It defines the required package sets:
  
  * CRAN: Shiny related packages plus `Seurat` and `ShinyCell`
* Bioconductor: `limma`, `edgeR`
If `do_markers = TRUE`, it also requires `presto`.
Behaviour:
  * If packages are missing and `install_missing = FALSE` it stops with a clear message listing missing packages.
* If `install_missing = TRUE` it installs missing CRAN packages with `install.packages()` and missing Bioconductor packages with `BiocManager::install()`.

3. Loads the Seurat object
You can provide either:
  
  * `seurat_obj` directly, or
* `seurat_rds` path to an `.rds` file
If both are missing it stops.

4. Creates the output directory
Ensures `out_dir` exists.

5. Validates and sets the RNA assay used for downstream steps
It checks `assay_rna` exists in the object and sets:
  
  * `DefaultAssay(seurat_obj) <- assay_rna`
This matters because variable features, marker detection, and exported matrices will follow this assay.

6. Optionally sets identities for clustering or grouping
If `ident_col` is provided, it sets:
  
  * `Idents(seurat_obj) <- meta.data[[ident_col]]`
This is important because many downstream Seurat functions rely on identities being correctly set.

7. Ensures assay keys exist
For each assay, it checks `Key(seurat_obj[[assay]])`.
If missing, it assigns a key like `rna_` or `atac_`.
This avoids collisions and prevents downstream plotting or feature naming issues.

8. Optionally computes variable features
If `do_variable_features = TRUE`, it runs:
  
  * `FindVariableFeatures(seurat_obj)`
This is mainly to ensure the object is in a sensible state for ShinyCell usage, especially if the object was created without variable features.

9. Optionally generates 3D UMAP reductions
If `do_umap3d = TRUE`, it will attempt to run 3D UMAP embeddings from one or more existing reductions (default is `pca`):
  
  For each reduction name in `umap3d_reductions`:
  
  * If the reduction exists in `seurat_obj@reductions`, it runs:
  
  * `RunUMAP(..., reduction = red, dims = umap3d_dims, n.components = 3, reduction.name = paste0(red, umap3d_name_suffix))`

Output effect:
  
  * Adds new reductions to the Seurat object such as `pca_umap3d` (if default suffix is used).
* These reductions can then be exported by ShinyCell and later used by ShinyCellPlus 3D tabs.

Important detail:
  
  * This does not “turn on” a 3D tab by itself. It only adds the embedding to the data so the 3D visualisation tab has something to plot later.

10. Optionally computes and writes marker lists to parquet
If `do_markers = TRUE` it builds a single marker table across one or more clustering resolution columns in `meta.data`.

How it decides which clusterings to use:
  
  * It finds metadata columns matching `markers_res_pattern` (default matches `res.`).

How markers are computed:
  
  * It extracts expression using `GetAssayData(..., layer = "data")` from the RNA assay.
* For each resolution column, it runs:
  
  * `presto::wilcoxauc(expr, clusters)`
* It appends a column `annotation` with the resolution name so multiple resolutions can live in one file.

Output behaviour:
  
  * Writes a parquet file to `markers_file` (default `out_dir/markergenes_lists.parquet`).
* If the parquet already exists and `markers_overwrite = FALSE`, it will skip and reuse the existing file.

Why this matters for the app:
  
  * The ShinyCellPlus marker tabs can load this parquet and display marker lists without recalculating them in Shiny.

11. Optionally exports raw counts to sc1counts.h5
If `do_counts_h5 = TRUE`, it writes a sparse raw counts matrix to an HDF5 file in compressed sparse column format (dgCMatrix CSC).

What it exports:
  
  * `counts <- GetAssayData(seurat_obj, assay = assay_rna, layer = counts_layer)` (default layer is `counts`)
* It enforces `dgCMatrix` and then exports the CSC slots:
  
  * `i`, `p`, `x`, `dims`
* plus `genes` and `cells`

Output behaviour:
  
  * Writes to `counts_h5_file` (default `out_dir/sc1counts.h5`)
* If file exists and `counts_overwrite = FALSE`, it skips and reuses it.
* It stores attributes describing the format, assay, and layer.

Why this matters for the app:
  
  * Tabs like pseudobulk can read counts on demand without loading a huge dense matrix into memory.

12. Creates a ShinyCell configuration and optionally builds the ShinyCell app folder
It always runs:
  
  * `scConf <- ShinyCell::createConfig(seurat_obj)`

If `do_make_app = TRUE`, it then runs:
  
  * `ShinyCell::makeShinyApp(seurat_obj, scConf, gene.mapping = gene_mapping, shiny.title = shiny_title, shiny.dir = out_dir)`

Output effect:
  
  * Writes the standard ShinyCell files into `out_dir` including:
  
  * `sc1conf.rds`, `sc1def.rds`, `sc1gene.rds`, `sc1meta.rds` (these are later loaded by `useShinyCellPlus()`).

13. Returns a structured summary invisibly
It returns a list containing:
  
  * the (possibly modified) `seurat_obj`
* `scConf`
* `out_dir`
* `markers_file` if markers were generated
* `counts_h5_file` if counts were generated

Key artefacts that future developers should know this function produces or enables:
  
  * Marker list file: `markergenes_lists.parquet` when `do_markers = TRUE`
* 3D UMAP reductions added to the Seurat object when `do_umap3d = TRUE`
* Counts H5 file: `sc1counts.h5` when `do_counts_h5 = TRUE`
* Standard ShinyCell output folder containing `sc1*.rds` when `do_make_app = TRUE`

If any of those are missing later in the app, it is because the corresponding `do_*` option was off, the file existed and overwrite was disabled, or the required reduction or metadata pattern was not found.


## useShinyCellPlus()

`useShinyCellPlus()` does not run analysis and does not compute markers or 3D UMAP. It generates and wires a runnable ShinyCellPlus application using already prepared data.

In plain terms, it takes:
  
  * `shiny.dir` → a folder containing prepared ShinyCell output files
* `shinycellplus.dir.src` → the ShinyCellPlus source modules
* configuration arguments (tabs, data type, title, etc.)

and produces a fully functional `app.R` that loads modules and launches the app.

Step by step, what it does:
  
  1. Validates and normalises paths
Ensures both the prepared data directory and the ShinyCellPlus source directory exist.

2. Requires an explicit app title
Stops execution if `app_title` is missing. Every app must be named.

3. Determines which tabs will be included
Tabs can be selected in two ways:
  
  • `data_type` selects predefined tab groups (RNA, RNA_ATAC, SPATIAL)
• `enabled_tabs` explicitly defines which tabs to include

If both are provided, it merges them.

The full set of available tab IDs that can be passed to `enabled_tabs` is:
  
  ```
c(
  "cellinfo_cellinfo",
  "cellinfo_geneexpr",
  "cellinfo3D_cellinfo3D",
  "cellinfo3D_geneexpr3D",
  "genecoex",
  "violin_boxplot",
  "proportions",
  "bubble_heatmap",
  "pseudobulk"
)
```

These correspond to modular components located in `modules/`.
Only tabs present in both `enabled_tabs` and the registered module list will be loaded.

4. Optionally disables legacy ui.R and server.R
If `disable_ui_server = TRUE`, existing `ui.R` and `server.R` files in `shiny.dir` are renamed to `.bak` so Shiny prioritises the generated `app.R`.

5. Copies the modular code
Copies the `modules/` folder from the ShinyCellPlus source directory into the target app directory.
If `overwrite_modules = TRUE`, it replaces any existing `modules/` folder.

6. Generates a new app.R
Writes an auto generated `app.R` file that:
  
  * Loads required libraries
* Defines shared theme and plotting utilities
* Loads prepared data files from `shiny.dir`:
  
  * `sc1conf.rds`
* `sc1def.rds`
* `sc1gene.rds`
* `sc1meta.rds`
* Detects whether `markergenes_lists.parquet` exists and passes its path to modules
* Sources all module files in `modules/`
* Registers and builds tabs dynamically
* Passes structured objects into each module server function:
  `sc1conf, sc1meta, sc1gene, sc1def, markers_list, assays, dir_inputs`

7. Optionally writes an rsconnect manifest
If `rsconnect.deploy = TRUE`, it writes deployment metadata but does not deploy automatically.

8. Returns the generated app path
Invisibly returns the path to the new `app.R`.

Important clarifications for developers:
  
  Marker lists
Markers are not computed here. If `markergenes_lists.parquet` exists in the prepared folder, its path is passed into modules for display and downstream use.

3D UMAP
3D visualisation depends on reductions prepared upstream. This function does not calculate embeddings. It only enables tabs such as:
  
  * `cellinfo3D_cellinfo3D`
* `cellinfo3D_geneexpr3D`
if they are requested and available.

sc1counts.h5
This function does not generate `sc1counts.h5`. It assumes expression data has already been prepared and stored upstream. Modules may access it using `dir_inputs`.

One sentence summary:
  `useShinyCellPlus()` assembles and configures a modular ShinyCellPlus app by copying source modules, generating `app.R`, loading prepared data objects, and activating the requested tab components.

## Each Module Structure

Each module is one R file inside ```modules```. each R script contains 4 blocks:
  
  - **functions**  (e.g., `scDRcell`)
- **UI function** (e.g., `scDRcell_ui`)
- **server function** (e.g., `scDRcell_server`)
- **registration** (e.g. `registration_tab`)

## Dependencies 

If you use additional libraries please add them to here or directly update the renv. Make sure you update the `useShinyCellPlus.R()` function

```
library(shiny)
library(shinyhelper)
library(shinyjs)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(hdf5r)
library(ggdendro)
library(gridExtra)
library(arrow)
library(rsconnect)
library(shinythemes)
library(shinydashboard)
library(tidyverse)
library(sortable)
library(plotly)
library(FlexDotPlot)
library(RColorBrewer)
library(ggforce)
library(limma) 
library(edgeR) 

```

## Create a new module

Modules in `modules` are ready to use. Modules under development live in `develop`.

Use the helper `createSCPModuleTemplate()` to create a new module file in `modules`. The helper function scans existing modules and avoids clashes with existing `register_tab()` ids and `sc1*` prefixes. Hence it is **crucial** to run `createSCPModuleTemplate()` in `modules`. If needed, you can move the generated file to `develop` afterwards.

```
createSCPModuleTemplate(
  module_dir  = "ShinyCellPlus/modules/",
  tab_id      = "scNewTab",
  module_name = "scNewTab",
  tab_title   = "New Analysis"
)
```

```
createSCPModuleTemplate <- function(
    module_dir,
    module_name,
    tab_id,
    tab_title,
    tab_header = tab_title,
    prefix_group = c("a","b","c","d","e"),
    overwrite = FALSE
)
```

The app loads modules automatically.

## Create a new module manually

If you want to create a module by hand, start by copying the relevant code from the original ShinyCell `ui.R` and `server.R`.

**1.** Find the tab you want to replicate
Search for an existing prefix, for example sc1a1, and copy the corresponding UI and server blocks for that tab.

**2.** In your new module file, keep the same structure
Include the following sections in order

- Function section
- UI section - create a UI function
- Server section - create a Server function
- Registration section

**3.** UI Function
At the top of the UI function, initialise the namespace

```
scDRnum_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
  ...
}
```
**3.1.** Namespace inputs in the UI
Wrap input ids with `ns()` in the UI. The ones that commonly get missed are `sc1a1`, `sc1a2`, `sc1b2`, etc.

**3.2.** Namespace ids inside the server where needed
Use `ns()` in places like `plotOutput()` and `checkboxGroupInput()` where you are declaring UI elements inside the server via `renderUI()`.

**3.3.** Update conditionalPanel() conditions
Use the namespaced input id inside sprintf, for example

```
condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2tog2"))
```

**3.4.** Make toggle labels user friendly
Rename toggle button text to something that reads clearly in the UI.

**4.** Server Function
At the top of the server function, initialise the session namespace and include the helper setup

```
scDRnum_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    ### For all tags and Server-side selectize 
    observe_helpers() 
    optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
    updateSelectizeInput(session, "sc1a1inp2", choices = names(sc1gene), server = TRUE, 
                         selected = sc1def$gene1, options = list( 
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1a3inp1", choices = names(sc1gene), server = TRUE, 
                         selected = sc1def$gene1, options = list( 
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1a3inp2", choices = names(sc1gene), server = TRUE, 
                         selected = sc1def$gene2, options = list( 
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1b2inp1", choices = names(sc1gene), server = TRUE, 
                         selected = sc1def$gene1, options = list( 
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1b2inp2", choices = names(sc1gene), server = TRUE, 
                         selected = sc1def$gene2, options = list( 
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1c1inp2", server = TRUE, 
                         choices = c(sc1conf[is.na(fID)]$UI,names(sc1gene)), 
                         selected = sc1conf[is.na(fID)]$UI[1], options = list( 
                           maxOptions = length(sc1conf[is.na(fID)]$UI) + 3, 
                           create = TRUE, persist = TRUE, render = I(optCrt))) 
    
    ...
  })
}
```

**4.1.** Note on server side selectize setup. If you are using server side selectize inputs, keep the observe_helpers() call and your updateSelectizeInput() setup near the start of the server function.

**5.** Register the tab
Your module must end with a register_tab() call like this

```
register_tab(
  id     = "cellinfo_geneexpr",
  title  = "CellInfo vs GeneExpr",
  ui     = scDRnum_ui,
  server = scDRnum_server
)

```

