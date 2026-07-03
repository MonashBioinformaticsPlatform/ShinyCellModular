# Prepare a Seurat object for ShinyCellModular

Takes a Seurat object from single cell experiments and prepares all
files needed to run a ShinyCellModular interactive Shiny app. Writes
HDF5 counts, ShinyCell config files, optional 3D UMAP, marker genes, and
motif data to `out_dir`.

## Usage

``` r
prepShinyCellModular(
  seurat_obj = NULL,
  seurat_rds = NULL,
  out_dir = "ShinyCellModular_app",
  shiny_title = "ShinyCellModular Intermediate",
  assays_selected = "RNA",
  ident_col = NULL,
  do_variable_features = TRUE,
  do_markers = FALSE,
  markers_file = NULL,
  markers_overwrite = FALSE,
  markers_res_pattern = "res\\.",
  do_umap3d = FALSE,
  umap3d_reductions = c("pca"),
  umap3d_dims = 1:30,
  umap3d_name_suffix = "_umap3d",
  do_counts_h5 = TRUE,
  counts_h5_file = NULL,
  counts_overwrite = TRUE,
  counts_layer = "counts",
  do_make_app = TRUE,
  gene_mapping = TRUE,
  install_missing = FALSE,
  verbose = TRUE,
  do_motifs = "auto",
  motifs_findmotifs = NULL,
  motifs_overwrite = TRUE,
  fragments_paths = NULL,
  custom_colors = NULL,
  default_genes = NULL,
  help = FALSE
)
```

## Arguments

- seurat_obj:

  Seurat object. Alternative to `seurat_rds`.

- seurat_rds:

  Path to a `.rds` Seurat object. Alternative to `seurat_obj`.

- out_dir:

  Output directory. Default: `'Files_ShinyCell'`.

- shiny_title:

  Title for the Shiny app.

- assays_selected:

  Assay(s) to process. e.g. `c('RNA', 'ATAC')`. Default: `'RNA'`.

- ident_col:

  Column to set as Idents. Default: `NULL` (uses existing).

- do_variable_features:

  Run `FindVariableFeatures` before creating the ShinyCell config.
  Default: `TRUE`.

- do_markers:

  Compute marker genes with presto. Default: `FALSE`.

- markers_file:

  Path for output markers parquet file. Default: auto.

- markers_overwrite:

  Overwrite existing markers file. Default: `FALSE`.

- markers_res_pattern:

  Regex pattern to find resolution columns. Default: `'res\\.'`.

- do_umap3d:

  Run 3D UMAP. Default: `FALSE`.

- umap3d_reductions:

  Reductions to use as input for 3D UMAP. Default: `c('pca')`.

- umap3d_dims:

  Dims to pass to UMAP (auto-capped to available). Default: `1:30`.

- umap3d_name_suffix:

  Suffix appended to the 3D UMAP reduction name. Default: `'_umap3d'`.

- do_counts_h5:

  Write raw counts to HDF5. Default: `TRUE`.

- counts_h5_file:

  Path for output H5 file. Default: auto.

- counts_overwrite:

  Overwrite existing H5 file. Default: `TRUE`.

- counts_layer:

  Seurat layer to use for counts. Default: `'counts'`.

- do_make_app:

  Run `makeShinyApp`. Default: `TRUE`.

- gene_mapping:

  Map gene names in ShinyCell. Default: `TRUE`.

- install_missing:

  Auto-install missing packages. Default: `FALSE`.

- verbose:

  Print progress messages. Default: `TRUE`.

- do_motifs:

  Extract motifs from ATAC assay. Runs automatically when ATAC is in
  `assays_selected` and the motifs slot is populated. Set to `FALSE` to
  skip. Default: `'auto'`.

- motifs_findmotifs:

  Output of `FindMotifs()`, adds enrichment scores. Default: `NULL`.

- motifs_overwrite:

  Overwrite existing motif files. Default: `TRUE`.

- fragments_paths:

  Optional named list of fragment file path overrides by index. e.g.
  `list('1' = '/path/to/sample1.tsv.gz')`. If `NULL`, paths are copied
  from the original paths in the object. Default: `NULL`.

- custom_colors:

  Named character vector of label -\> hex color to override ShinyCell
  default colors. Extra names not in the data are ignored. Default:
  `NULL`.

- default_genes:

  Character vector of gene names to set as defaults in the app (e.g.
  `c('CD4', 'CD8A')`). If `NULL`, ShinyCell picks its own defaults.
  Default: `NULL`.

- help:

  Print the help message and return. Default: `FALSE`.

## Value

Invisibly returns `NULL`. Writes output files to `out_dir`.

## Examples

``` r
if (FALSE) { # \dontrun{
prepShinyCellModular(
  seurat_rds = "seurat_object.rds",
  out_dir    = "my_app_files",
  do_umap3d  = TRUE,
  do_markers = TRUE
)
} # }
```
