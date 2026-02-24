# ShinyCellPlus

ShinyCellPlus is a modular, enhanced version of ShinyCell developed at the Monash Genomics and Bioinformatics Platform (MGBP). It supports large scRNAseq and multimodal datasets with fast on-demand HDF5 access, extended visualisations, improved filtering, and publication-ready plots. Its modular structure makes it flexible, scalable, and easy to customise.

## Features

- Modular UI and server structure
- Supports scRNAseq, ATAC, and multimodal datasets
- Fast HDF5 on-demand loading
- Publication‑ready plots (PNG/PDF export)
- Extended visualisation tabs (UMAP, violin, bubble, heatmap, coexpression, ROC)
- Dynamic selectize-based gene input
- Cell subsetting and conditional plotting
- Easy integration with new modules via a registry system

## Structure

Each module is one R file containing:

- A plotting **function** (e.g., `scDRcell`)
- A **UI function** (e.g., `scDRcellgene_ui`)
- A **server function** (e.g., `scDRcellgene_server`)

<<<<<<< HEAD

## Depencecies

```
#load the Project 
renv::restore(prompt = FALSE)
renv::snapshot(prompt = FALSE)
renv::status()

renv::restore()

# Make sure you have all this packages
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
library(plotly)

```

```
code/modules/
<<<<<<< HEAD

scDRnum.R <- "CellInfo vs GeneExpr" - 1 Tab


bubbleplot replacement:

https://github.com/Simon-Leonard/FlexDotPlot

=======
>>>>>>> 4cf02dd (cellinfo_geneexp working and README.md added)
```

A central registry lists all modules:

```
code/registry.R
```

The app loads modules automatically.

## Usage

Place your ShinyCell-processed files:

```
sc1conf.rds
sc1meta.rds
sc1gene.rds
sc1def.rds
sc1gexpr.h5
```

Then run:

```r
<<<<<<< HEAD
source("../useShinyCellPlus.R")
useShinyCellPlus(
    shiny.dir="Files_ShinyCell/",
    shinycellplus.dir.src="~/Dropbox/BioPlatform/ShinyCellPlus_devel/ShinyCellPlus/",
    rsconnect.deploy = FALSE,
    enabled_tabs = "cellinfo_geneexpr",
    overwrite_code = TRUE
)

This step created a app.R file. This file is the only file you need to run to run your ShinyCellPlus app.

=======
shiny::runApp("app.R")
>>>>>>> 4cf02dd (cellinfo_geneexp working and README.md added)
```


# TEMPLATE for new tabs

to create a new module copy paste all related functions from the server.r. 
Copy paste the section of ui.r and server.r that correspond to that tab . 

in the same file make a function section, UI section, server Section, and Registration Section
create a ui and server function. 



add ns() to all inputs in ui. namespaced. search "sc1a1/a2/b2" those are usually the ones that need ns()
in server add ns() to   plotOutput and checkboxGroupInput


change condition panel like:
```
condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2tog2"))
```

change the work Toggle  for more user frienly language

this at the start of UI:

```
scDRnum_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
```


add this to the start of server:
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
    
```

############################################### Registration #################################################
```
register_tab(
  id     = "cellinfo_geneexpr",
  title  = "CellInfo vs GeneExpr",
  ui     = scDRnum_ui,
  server = scDRnum_server
)
```





## Adding a Module

<<<<<<< HEAD
1. Create the app source under `code/modules/`
=======
1. Create a file under `code/modules/`
>>>>>>> 4cf02dd (cellinfo_geneexp working and README.md added)
2. Define `module_ui`, `module_server`, and any helper functions
3. Add an entry to `tab_registry` in `code/registry.R`
4. Include the module key in `enabled_tabs` in `app.R`

## Credits

- Built on **ShinyCell** (https://github.com/SGDDNB/ShinyCell)
- MGBP extended version **ShinyCellPlus** (https://github.com/MonashBioinformaticsPlatform/ShinyCellPlus)

## License

MIT

