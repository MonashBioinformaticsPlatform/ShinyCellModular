# prepShinyCellModular()

`prepShinyCellModular()` is the preparation step that builds everything the ShinyCellModular app will need. It takes a Seurat object, processes one or more assays (RNA, ATAC), and produces a ready-to-use directory on disk containing `.rds`/`.h5`/`.parquet` files that `useShinyCellModular()` and the modular tabs consume.

Each assay in `assays_selected` gets its own output subfolder (`out_dir/RNA/`, `out_dir/ATAC/`), so a single call can prepare a multiome object end to end.

---

## Function signature

```r
prepShinyCellModular(
    seurat_obj            = NULL,
    seurat_rds            = NULL,
    out_dir               = "ShinyCellModular_app",
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
    custom_colors         = NULL,
    default_genes         = NULL,
    help                  = FALSE
)
```

Call `prepShinyCellModular(help = TRUE)` at any time to print the same argument summary from R (also printed automatically if neither `seurat_obj` nor `seurat_rds` is supplied).

---

## Arguments

| Argument | Default | Description |
|---|---|---|
| `seurat_obj` | `NULL` | Seurat object. Alternative to `seurat_rds`. |
| `seurat_rds` | `NULL` | Path to a `.rds` Seurat object. Alternative to `seurat_obj`. |
| `out_dir` | `"ShinyCellModular_app"` | Output directory. |
| `shiny_title` | `"ShinyCellModular Intermediate"` | Title for the Shiny app. |
| `assays_selected` | `"RNA"` | Assay(s) to process, e.g. `c("RNA","ATAC")`. |
| `ident_col` | `NULL` | Column to set as `Idents`. `NULL` keeps the existing idents. |
| `do_variable_features` | `TRUE` | Run `FindVariableFeatures()` before `createConfig()`. |
| `do_markers` | `FALSE` | Compute marker genes with `presto::wilcoxauc`. RNA only. |
| `markers_file` | `NULL` | Output markers parquet path. Auto: `out_dir/RNA/markergenes_lists.parquet`. |
| `markers_overwrite` | `FALSE` | Overwrite an existing markers file. |
| `markers_res_pattern` | `"res\\."` | Regex used to find clustering-resolution columns in `meta.data`. |
| `do_umap3d` | `FALSE` | Run a 3D UMAP (`n.components = 3`) for each reduction in `umap3d_reductions`. |
| `umap3d_reductions` | `c("pca")` | Reductions used as input for the 3D UMAP. |
| `umap3d_dims` | `1:30` | Dims passed to `RunUMAP()`. |
| `umap3d_name_suffix` | `"_umap3d"` | Suffix appended to the new reduction name, e.g. `pca_umap3d`. |
| `do_counts_h5` | `TRUE` | Write raw/sparse counts to an HDF5 file (`sc1counts.h5`). |
| `counts_h5_file` | `NULL` | Output H5 path. Auto: `out_dir/RNA/sc1counts.h5`. |
| `counts_overwrite` | `TRUE` | Overwrite an existing H5 file. |
| `counts_layer` | `"counts"` | Seurat layer read via `GetAssayData()`. |
| `do_make_app` | `TRUE` | Run `ShinyCell::createConfig()` + `ShinyCell::makeShinyApp()`. |
| `gene_mapping` | `TRUE` | Passed through to `makeShinyApp(gene.mapping = ...)`. |
| `install_missing` | `FALSE` | Auto-install missing CRAN/Bioconductor/GitHub dependencies instead of stopping. |
| `verbose` | `TRUE` | Print progress messages via `message()`. |
| `do_motifs` | `"auto"` | Extract motifs from the ATAC assay. `"auto"` runs only if a motifs slot is populated; `TRUE`/`FALSE` force it on/off. |
| `motifs_findmotifs` | `NULL` | Output of `Signac::FindMotifs()`; joined into motif metadata as enrichment scores. |
| `motifs_overwrite` | `TRUE` | Overwrite existing `sc1motifs.rds` / `sc1motifs_meta.parquet`. |
| `fragments_paths` | `NULL` | Named list of fragment-file overrides keyed by fragment index, e.g. `list("1" = "/path/sample1.tsv.gz")`. Used when the original path stored in the object is no longer valid. |
| `custom_colors` | `NULL` | Named character vector (`level = "#hexcolor"`) patched into `sc1conf.rds` after `makeShinyApp()` writes it. Unmatched names are ignored. |
| `default_genes` | `NULL` | Character vector of genes to preselect in the app (`default.multigene`). |
| `help` | `FALSE` | Print the help message and return, without processing anything. |

**Return value:** invisibly returns a list — `seurat_obj`, `out_dir`, `markers_file` (if `do_markers`), `counts_h5_file` (if `do_counts_h5`). Files are the real output; the return value is a convenience for chaining in scripts.

---

## Step by step

### 1. Dependency check
Builds a list of required CRAN, Bioconductor, and GitHub-only packages (`ShinyCell`, `ggseqlogo`, `FlexDotPlot`, and `presto` if `do_markers = TRUE`; `Signac` is added automatically if `"ATAC" %in% assays_selected`). Missing packages stop execution with an install hint unless `install_missing = TRUE`, in which case they're installed automatically (`install.packages()`, `BiocManager::install()`, `remotes::install_github()`).

### 2. Load the Seurat object
Exactly one of `seurat_obj` / `seurat_rds` must resolve to a Seurat object. If both are `NULL` (and `help` wasn't explicitly set), `help` is forced to `TRUE` and the function prints `helpMessage` and returns. Passing a file path to `seurat_obj` or an in-memory object to `seurat_rds` triggers a `warning()` pointing at the likely mistake.

### 3. Validate `assays_selected`
Checks every entry exists in `Seurat::Assays(seurat_obj)`; stops with the requested vs. available assay names otherwise.

### 4. Ensure assay keys exist
For any assay with a missing/empty `Key()`, sets one automatically (`tolower(assay) + "_"`) so downstream Seurat/Signac calls don't fail.

### 5. Main loop over `assays_selected`
For each `active_assay` (`DefaultAssay(seurat_obj) <- active_assay`):

**RNA branch** (`active_assay == "RNA"`):
- Creates `out_dir/RNA/`.
- Optionally sets `Idents()` from `ident_col`.
- Optionally adds 3D UMAP reductions (`do_umap3d`), skipping any reduction not present with a message.
- Optionally computes markers per resolution column matched by `markers_res_pattern` using `presto::wilcoxauc`, skipping resolutions with fewer than 2 groups, and writes the combined table to `markers_file` via `arrow::write_parquet()`.

**ATAC branch** (`active_assay == "ATAC"`):
- Creates `out_dir/ATAC/`.
- Resolves whether to run motif extraction from `do_motifs` (`"auto"` checks the `@motifs` slot; otherwise `isTRUE(do_motifs)`), and if so extracts PWMs + names into `sc1motifs.rds` and `sc1motifs_meta.parquet`, optionally joined with `motifs_findmotifs` enrichment scores.
- Always extracts static ATAC objects: `sc1annotation.rds`, `sc1peaks.rds`, `sc1links.rds` (each skipped individually with a message if the slot is empty), and copies fragment files into `out_dir/ATAC/fragments/index<N>/`, using `fragments_paths` overrides for any file that's missing from the object's original path. Writes `sc1fragmentpaths.rds`.

**All assays:**
- If `do_make_app = TRUE`: optionally runs `FindVariableFeatures()`, then `ShinyCell::createConfig()` and `ShinyCell::makeShinyApp()` into `out_dir/<assay>/`. If `custom_colors` is supplied, patches matching factor-level colours into the `sc1conf.rds` that `makeShinyApp()` just wrote (since `makeShinyApp()` ignores any pre-patched config passed in).
- If `do_counts_h5 = TRUE`: writes the sparse counts matrix (CSC format: `i`, `p`, `x`, `dims`, `genes`, `cells`) to `out_dir/<assay>/sc1counts.h5` via `hdf5r`.

---

## Notes

- `custom_colors` passed via `source("colors.R")` returns a list, not a named vector — the function detects this, extracts `$value` automatically, and warns you to pass `custom_colors = source("colors.R")$value` directly next time.
- Fragment files that can't be found anywhere (neither the original path nor a `fragments_paths` override) are reported with a ready-to-paste `fragments_paths = list(...)` snippet in the message.
- `install_missing = TRUE` also runs `setRepositories(ind = 1:3)` before `install.packages()` so Bioconductor dependencies (e.g. for `Signac`) resolve correctly.
