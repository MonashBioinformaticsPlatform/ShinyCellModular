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

Modules are stored in:

```
code/modules/
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
shiny::runApp("app.R")
```

## Adding a Module

1. Create a file under `code/modules/`
2. Define `module_ui`, `module_server`, and any helper functions
3. Add an entry to `tab_registry` in `code/registry.R`
4. Include the module key in `enabled_tabs` in `app.R`

## Credits

- Built on **ShinyCell** (https://github.com/SGDDNB/ShinyCell)
- MGBP extended version **ShinyCellPlus** (https://github.com/MonashBioinformaticsPlatform/ShinyCellPlus)

## License

MIT

