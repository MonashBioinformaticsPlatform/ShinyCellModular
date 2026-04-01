# prepShinyCellModular()

`prepShinyCellModular()` is the preparation step that builds everything the ShinyCellModular app will need. It takes a Seurat object, processes one or more assays (RNA, ATAC), and produces a ready-to-use directory on disk containing `.rds` files and optional extras that `useShinyCellModular()` and the modular tabs can consume.

---

## Function signature

```r
prepShinyCellModular(
    seurat_obj            = NULL,
    seurat_rds            = NULL,
    out_dir               = "Files_ShinyCell",
    shiny_title           = "ShinyCellModular Intermediate",
    assays_selected       = "RNA",
    ident_col             = NULL,
    do_variable_features  = TRUE,
    do_markers            = FALSE,
    markers_file          = NULL,
    markers_overwrite     = FALSE,
    markers_res_pattern   = "res\\.",
    do_umap3d             = FALSE,
    umap3d_reductions     = c("pca"),
    umap3d_dims           = 1:30,
    umap3d_name_suffix    = "_umap3d",
    do_counts_h5          = TRUE,
    counts_h5_file        = NULL,
    counts_overwrite      = TRUE,
    counts_layer          = "counts",
    do_make_app           = TRUE,
    gene_mapping          = TRUE,
    install_missing       = FALSE,
    verbose               = TRUE,
    do_motifs             = "auto",
    motifs_findmotifs     = NULL,
    motifs_overwrite      = TRUE,
    fragments_paths       = NULL,
    fragments_sample_col  = NULL
)
```

---

## Step by step

### 1. Sets file path defaults
If not provided:
- `markers_file` defaults to `out_dir/markergenes_lists.parquet`
- `counts_h5_file` defaults to `out_dir/sc1counts.h5`

### 2. Checks dependencies and optionally installs them
Defines required package sets:
- CRAN: Shiny-related packages, `Seurat`, `ShinyCell`
- Bioconductor: `limma`, `edgeR`
- If `do_markers = TRUE`: also requires `presto`

Behaviour:
- If packages are missing and `install_missing = FALSE`, stops with a clear message listing missing packages
- If `install_missing = TRUE`, installs missing CRAN packages with `install.packages()` and Bioconductor packages with `BiocManager::install()`

### 3. Loads the Seurat object
Either from `seurat_obj` directly or from `seurat_rds` path. Stops if neither is provided.

### 4. Creates the output directory
Ensures `out_dir` exists.

### 5. Validates assays_selected
Checks that every assay in `assays_selected` exists in the Seurat object. Stops with a clear message listing missing assays if any are not found.

### 6. Ensures assay keys exist
For every assay in the object, checks `Key()`. If missing, assigns a key like `rna_` or `atac_`. This prevents downstream collisions in ShinyCell and plotting functions.

### 7. Loops over each assay in assays_selected
All subsequent steps run inside this loop. `DefaultAssay()` is set at the top of each iteration.

---

## RNA-specific steps

Only run when the active assay is `RNA`.

**Optionally sets identities**
If `ident_col` is provided, sets `Idents(seurat_obj) <- meta.data[[ident_col]]`. Many downstream Seurat functions rely on identities being correctly set.

**Optionally computes 3D UMAP**
If `do_umap3d = TRUE`, runs `RunUMAP()` with `n.components = 3` for each reduction listed in `umap3d_reductions`. Results are stored as `<reduction><umap3d_name_suffix>` (default: `pca_umap3d`). These reductions are exported by ShinyCell and consumed by the 3D visualisation tabs. This does not activate a 3D tab by itself — it only adds the embedding to the data.

**Optionally computes marker genes**
If `do_markers = TRUE`:
- Finds metadata columns matching `markers_res_pattern` (default matches `res.`)
- Extracts expression with `GetAssayData(..., layer = "data")`
- Runs `presto::wilcoxauc()` per resolution column
- Appends an `annotation` column with the resolution name so multiple resolutions live in one file
- Writes a single parquet to `markers_file`
- Skips if file exists and `markers_overwrite = FALSE`

---

## ATAC-specific steps

Only run when the active assay is `ATAC`. Creates a subfolder `out_dir/ATAC/`.

**Motif extraction**
Controlled by `do_motifs`:
- `"auto"` (default): runs if a motifs slot is detected in the assay
- `TRUE`: always runs
- `FALSE`: always skips

When run:
- Saves PWM list to `sc1motifs.rds`
- Saves motif metadata to `sc1motifs_meta.parquet`
- If `motifs_findmotifs` is provided (output of `FindMotifs()`), enrichment scores are joined into the metadata
- Skips if files exist and `motifs_overwrite = FALSE`

**Static ATAC object extraction**
Always runs for ATAC. Saves the following to `out_dir/ATAC/`:

| File | Contents | Skipped if |
|---|---|---|
| `sc1annotation.rds` | Genome annotation | No annotation slot found |
| `sc1peaks.rds` | Peak ranges | No ranges slot found |
| `sc1links.rds` | Peak-to-gene links | No links found |
| `sc1fragmentpaths.rds` | Fragment file paths and cell barcodes | No fragment files found |

Fragment files are copied into `out_dir/ATAC/fragments/`. If original paths are missing, a clear message is shown with the `fragments_paths` override syntax:

```r
fragments_paths = list(
  "1" = "/path/to/sample1.tsv.gz",
  "2" = "/path/to/sample2.tsv.gz"
)
```

---

## Steps that run for all assays

**Optionally creates ShinyCell files**
If `do_make_app = TRUE`:
- Runs `FindVariableFeatures()` if `do_variable_features = TRUE`
- Runs `ShinyCell::createConfig()` and `ShinyCell::makeShinyApp()`
- Output goes to `out_dir` for RNA, or `out_dir/<ASSAY>/` for other assays
- Produces `sc1conf.rds`, `sc1def.rds`, `sc1gene.rds`, `sc1meta.rds`

**Optionally exports raw counts to H5**
If `do_counts_h5 = TRUE`:
- Extracts `GetAssayData(..., layer = counts_layer)` (default: `counts`)
- Enforces `dgCMatrix` and writes CSC slots (`i`, `p`, `x`, `dims`) plus `genes` and `cells`
- Output goes to `counts_h5_file` for RNA, or `out_dir/<ASSAY>/sc1counts.h5` for other assays
- Skips if file exists and `counts_overwrite = FALSE`
- Tabs like `pseudobulk` use this to read counts on demand without loading a dense matrix

---

## Return value

Returns invisibly:

```r
list(
    seurat_obj     = seurat_obj,
    out_dir        = out_dir,
    markers_file   = if (do_markers)   markers_file   else NULL,
    counts_h5_file = if (do_counts_h5) counts_h5_file else NULL
)
```

---

## Key artefacts produced

| Artefact | Location | Requires |
|---|---|---|
| `sc1conf.rds`, `sc1def.rds`, `sc1gene.rds`, `sc1meta.rds` | `out_dir/` | `do_make_app = TRUE` |
| `markergenes_lists.parquet` | `out_dir/` | `do_markers = TRUE` |
| `sc1counts.h5` | `out_dir/` | `do_counts_h5 = TRUE` |
| `pca_umap3d` reduction in Seurat object | (in memory, exported via ShinyCell) | `do_umap3d = TRUE` |
| `sc1conf.rds` etc. for ATAC | `out_dir/ATAC/` | `"ATAC" %in% assays_selected`, `do_make_app = TRUE` |
| `sc1counts.h5` for ATAC | `out_dir/ATAC/` | `"ATAC" %in% assays_selected`, `do_counts_h5 = TRUE` |
| `sc1motifs.rds`, `sc1motifs_meta.parquet` | `out_dir/ATAC/` | `"ATAC" %in% assays_selected`, motifs slot present |
| `sc1annotation.rds`, `sc1peaks.rds`, `sc1links.rds` | `out_dir/ATAC/` | `"ATAC" %in% assays_selected`, slots present |
| `sc1fragmentpaths.rds` + copied fragment files | `out_dir/ATAC/fragments/` | `"ATAC" %in% assays_selected`, fragments present |

If any artefact is missing later in the app, it is because the corresponding `do_*` option was off, the file existed and overwrite was disabled, the required slot was not found in the object, or the reduction or metadata pattern was not matched.


***

# useShinyCellModular()

`useShinyCellModular()` does not run analysis and does not compute markers, embeddings, or counts. It generates and wires a runnable ShinyCellModular application using already prepared data.

In plain terms, it takes:
- `shiny.dir` — a folder containing prepared ShinyCell output files
- `shinycellmodular.dir.src` — the ShinyCellModular source directory containing `modules/`
- configuration arguments (tabs, data type, title, etc.)

and produces a fully functional `app.R` that loads modules dynamically and launches the app.

---

## Function signature

```r
useShinyCellModular(
    shiny.dir,
    shinycellmodular.dir.src,
    rsconnect.deploy  = FALSE,
    data_type         = c("RNA", "RNA_ATAC", "SPATIAL"),
    enabled_tabs      = NULL,
    overwrite_modules = FALSE,
    disable_ui_server = TRUE,
    app_title         = NULL
)
```

---

## Step by step

### 1. Validates and normalises paths
Ensures both `shiny.dir` and `shinycellmodular.dir.src` exist. Stops if either is missing.

### 2. Requires an explicit app title
Stops execution if `app_title` is missing or `NULL`. Every app must be named.

### 3. Determines which tabs will be included
The available tabs are discovered dynamically by scanning the `modules/<data_type>/` subfolder inside `shinycellmodular.dir.src`. This means the tab catalogue is driven by what module files are actually present on disk, not a hardcoded list.

Tabs can be selected in two ways:
- `data_type` selects all tabs found in the corresponding subfolder
- `enabled_tabs` explicitly defines which subset of those tabs to include

Rules:
- If only `data_type` is provided, all tabs for that type are included
- If only `enabled_tabs` is provided, `data_type` defaults to `"RNA"` with a warning
- If both are provided, `enabled_tabs` is validated against the tabs available for `data_type`. Any tab not belonging to that type causes an error with a clear message listing valid options
- If neither is provided, execution stops

### 4. Optionally disables legacy ui.R and server.R
If `disable_ui_server = TRUE` and `ui.R` or `server.R` exist in `shiny.dir`, they are renamed to `.bak`. This ensures Shiny picks up the generated `app.R` rather than the legacy files. A warning is shown when this happens.

### 5. Copies the requested module files
Copies only the module files corresponding to `enabled_tabs` from `shinycellmodular.dir.src/modules/<data_type>/` into `shiny.dir/modules/`. Only the requested tabs are copied — not the entire modules folder — to avoid bloat as the project grows.

If `overwrite_modules = TRUE`, the existing `modules/` folder in `shiny.dir` is removed and replaced. A warning is shown. Any local modifications inside `modules/` will be lost.

### 6. Generates app.R
Writes an auto-generated `app.R` into `shiny.dir` that:
- Loads all required libraries
- Defines shared theme and plotting utilities (`cList`, `pList`, `sList`, `sctheme`, `g_legend`)
- Loads prepared data files from `shiny.dir`:
  - `sc1conf.rds`, `sc1def.rds`, `sc1gene.rds`, `sc1meta.rds`
- Loads ATAC files from `shiny.dir/ATAC/` if that folder exists:
  - `sc1conf.rds`, `sc1def.rds`, `sc1gene.rds`, `sc1meta.rds` (ATAC versions)
  - `sc1fragmentpaths.rds`, `sc1annotation.rds`, `sc1peaks.rds`, `sc1links.rds`
  - All loaded with `tryCatch` — missing files return `NULL` gracefully
- Detects whether `markergenes_lists.parquet` exists and passes its path to modules as `markers_list`, or `NULL` if absent
- Sources all module files in `modules/` dynamically using `source()`
- Registers tabs via `register_tab()` calls inside each module
- Builds the UI as a `navbarPage` with one `tabPanel` per enabled tab
- Wires the server by calling each tab's server function with only the arguments it declares (matched via `formals()`)

Arguments passed to each module server (if declared in its signature):

```r
id, sc1conf, sc1meta, sc1gene, sc1def,
sc1conf_atac, sc1def_atac, sc1gene_atac, sc1meta_atac,
sc1fragmentpaths, sc1annotation, sc1peaks, sc1links,
markers_list, assays, dir_inputs
```

### 7. Optionally writes an rsconnect manifest
If `rsconnect.deploy = TRUE`:
- Runs `rsconnect::writeManifest()` in `shiny.dir`
- Patches the manifest to fix paths for `.rds`, `.h5`, and `.parquet` files
- Does not deploy automatically — the manifest is written for manual deployment

### 8. Returns the app path
Returns the path to the generated `app.R` invisibly.

---

## Tab discovery

 `useShinyCellModular()` discovers available tabs by scanning the `modules/` directory structure:

```
shinycellmodular.dir.src/
  modules/
    RNA/
      bubble_heatmap.R
      violin_boxplot.R
      ...
    RNA_ATAC/
      ...
    SPATIAL/
      ...
```

Each subfolder corresponds to a `data_type`. The tab IDs are the filenames without extension. This means adding a new tab only requires dropping a new module file into the right subfolder — no changes to `useShinyCellModular()` are needed.

---

## Important clarifications

**Marker lists**
Markers are not computed here. If `markergenes_lists.parquet` exists in `shiny.dir`, its path is passed to modules. If absent, `markers_list` is `NULL` and tabs must handle that gracefully.

**3D UMAP**
3D tabs depend on reductions prepared by `prepShinyCellModular()` with `do_umap3d = TRUE`. This function does not compute embeddings. It only enables 3D tabs if they are present in the module folder and requested via `enabled_tabs`.

**ATAC data**
ATAC-specific objects (`sc1fragmentpaths`, `sc1annotation`, `sc1peaks`, `sc1links`) are loaded from `shiny.dir/ATAC/` if that folder exists. All loads use `tryCatch` so missing files return `NULL` without crashing the app. Tabs that require these objects must check for `NULL` themselves.

**sc1counts.h5**
Not loaded globally. Modules access it on demand via `dir_inputs`, which points to `shiny.dir`.

**Module argument matching**
The server wiring uses `formals(srv)` to match arguments. A module server will only receive the arguments it explicitly declares. This means modules are forward-compatible — new arguments can be added to the wiring without breaking existing modules that do not declare them.

---

## One sentence summary

`useShinyCellModular()` assembles a modular ShinyCellModular app by discovering available tabs from the modules folder, copying requested modules, generating `app.R` with dynamic tab loading and data wiring, and optionally preparing an rsconnect manifest.


# Considerations to pass `enabled_tabs`
#### Tab Catalogue

Available tabs for `useShinyCellModular()`. Tab IDs are passed to `enabled_tabs`. Tabs are discovered dynamically from the `modules/<data_type>/` subfolder — the tab ID matches the module filename without extension.

## RNA

| Tab id (`enabled_tabs`)   | Tab title (UI)           | Module file                | What it contains                                                                                  | What you need in `prepShinyCellModular`                                                                                         | Developed |
| ------------------------- | ------------------------ | -------------------------- | ------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- | --------- |
| `cellinfo_cellinfo`       | CellInfo vs CellInfo     | `cellinfo_cellinfo.R`      | 2D embedding coloured by metadata, with optional grouping or splitting by a second metadata field | Standard ShinyCell outputs (`sc1conf.rds`, `sc1meta.rds`, `sc1gene.rds`, `sc1def.rds`). Metadata columns must exist before prep | Yes       |
| `cellinfo_geneexpr`       | CellInfo vs GeneExpr     | `cellinfo_geneexpr.R`      | 2D embedding with gene expression overlay (feature style plots)                                   | Standard ShinyCell outputs. `gene_mapping = TRUE` recommended if you rely on gene aliases                                       | Yes       |
| `cellinfo3D_cellinfo3D`   | CellInfo3D vs CellInfo3D | `cellinfo3D_cellinfo3D.R`  | Interactive 3D embedding coloured by metadata                                                     | `do_umap3d = TRUE` so 3D reductions exist (default from PCA). Ensure `umap3d_reductions` exists in the object                  | Yes       |
| `cellinfo3D_geneexpr3D`   | CellInfo3D vs GeneExpr3D | `cellinfo3D_geneexpr3D.R`  | Interactive 3D embedding with gene expression overlay                                             | `do_umap3d = TRUE` so 3D reductions exist                                                                                      | Yes       |
| `genecoex`                | Gene Coexpression        | `genecoex.R`               | Coexpression visualisation for selected genes across cells or groups                              | Standard ShinyCell outputs                                                                                                      | Yes       |
| `violin_boxplot`          | Violin / BoxPlot         | `violin_boxplot.R`         | Violin and boxplots for gene expression or metadata across groups                                 | Standard ShinyCell outputs                                                                                                      | Yes       |
| `proportions`             | Cell Proportions         | `proportions.R`            | Composition summaries across groups (for example cluster proportions per sample)                  | Standard ShinyCell outputs. Requires grouping metadata such as sample and cluster                                                | Yes       |
| `bubble_heatmap`          | Bubble Plot / Heatmap    | `bubble_heatmap.R`         | Bubble plot and heatmap summaries for gene sets across groups with size and colour encodings       | Standard ShinyCell outputs                                                                                                      | Yes       |
| `pseudobulk`              | Pseudobulk DE            | `pseudobulk.R`             | Pseudobulk aggregation and DE workflow based on raw counts                                        | `do_counts_h5 = TRUE` to create `sc1counts.h5`. Requires `edgeR` and `limma`                                                   | Yes       |

---

## RNA_ATAC

Multimodal tabs for datasets with both RNA and ATAC assays. These are separate module files from the RNA equivalents and carry the `_multi` suffix. Both assays must be prepared with `prepShinyCellModular` using `assays_selected = c("RNA", "ATAC")`.

| Tab id (`enabled_tabs`)         | Tab title (UI)           | Module file                       | What it contains                                                                                  | What you need in `prepShinyCellModular`                                                                                                        | Developed |
| ------------------------------- | ------------------------ | --------------------------------- | ------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- | --------- |
| `cellinfo_cellinfo_multi`       | CellInfo vs CellInfo     | `cellinfo_cellinfo_multi.R`       | 2D embedding coloured by metadata, multimodal version                                             | Standard ShinyCell outputs for RNA and ATAC (`out_dir/ATAC/`)                                                                                  | Yes       |
| `cellinfo_geneexpr_multi`       | CellInfo vs GeneExpr     | `cellinfo_geneexpr_multi.R`       | 2D embedding with gene expression overlay, multimodal version                                     | Standard ShinyCell outputs for RNA and ATAC                                                                                                    | Yes       |
| `cellinfo3D_cellinfo3D_multi`   | CellInfo3D vs CellInfo3D | `cellinfo3D_cellinfo3D_multi.R`   | Interactive 3D embedding coloured by metadata, multimodal version                                 | `do_umap3d = TRUE`. Standard ShinyCell outputs for RNA and ATAC                                                                                | Yes       |
| `cellinfo3D_geneexpr3D_multi`   | CellInfo3D vs GeneExpr3D | `cellinfo3D_geneexpr3D_multi.R`   | Interactive 3D embedding with gene expression overlay, multimodal version                         | `do_umap3d = TRUE`. Standard ShinyCell outputs for RNA and ATAC                                                                                | Yes       |
| `genecoex_multi`                | Gene Coexpression        | `genecoex_multi.R`                | Coexpression visualisation for selected genes, multimodal version                                 | Standard ShinyCell outputs for RNA and ATAC                                                                                                    | Yes       |
| `violin_boxplot_multi`          | Violin / BoxPlot         | `violin_boxplot_multi.R`          | Violin and boxplots for gene expression or metadata, multimodal version                           | Standard ShinyCell outputs for RNA and ATAC                                                                                                    | Yes       |
| `proportions_multi`             | Cell Proportions         | `proportions_multi.R`             | Composition summaries across groups, multimodal version                                           | Standard ShinyCell outputs for RNA and ATAC. Requires grouping metadata such as sample and cluster                                             | Yes       |
| `bubble_heatmap_multi`          | Bubble Plot / Heatmap    | `bubble_heatmap_multi.R`          | Bubble plot and heatmap summaries, multimodal version                                             | Standard ShinyCell outputs for RNA and ATAC                                                                                                    | Yes       |
| `pseudobulk_multi`              | Pseudobulk DE            | `pseudobulk_multi.R`              | Pseudobulk aggregation and DE workflow, multimodal version                                        | `do_counts_h5 = TRUE` for RNA and ATAC. Requires `edgeR` and `limma`                                                                          | Yes       |
| `coverage_plot`                 | Coverage Plot            | `coverage_plot.R`                 | ATAC signal coverage tracks over genomic regions                                                  | `"ATAC" %in% assays_selected`. Fragment files accessible via `sc1fragmentpaths.rds`. Annotation (`sc1annotation.rds`) recommended             | Yes       |
| `motif_plot`                    | Motif Plot               | `motif_plot.R`                    | Transcription factor motif activity visualisation                                                 | `"ATAC" %in% assays_selected`. Motifs slot present in ATAC assay so `prepShinyCellModular` can extract `sc1motifs.rds`                         | Yes       |

---

## CROPSEQ

No modules developed yet.

---

## SPATIAL

No modules developed yet.
