
# Functions Details

```
prepShinyCellPlus <- function(
    seurat_obj = NULL,
    seurat_rds = NULL,
    out_dir = "Files_ShinyCell",
    shiny_title = "ShinyCellPlus Intermediate",
    assay_rna = "RNA",
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
    verbose = TRUE
) 
```

***

```
useShinyCellPlus <- function(
    shiny.dir, # files from shinycell are
    shinycellplus.dir.src, # modules where shinycellplus 
    rsconnect.deploy = FALSE, # do you want to publish in rsconnect
    data_type = c("RNA", "RNA_ATAC", "SPATIAL"), # what predetermine tabs you want
    enabled_tabs = NULL, # what tabs you want
    overwrite_modules = FALSE, # overwrite modules
    disable_ui_server = TRUE, # this disables the existing ui.R and server.r
    app_title=NULL
) 

```


### Considerations to pass `enabled_tabs`

  | Tab id (enabled_tabs)   | Tab title (UI)           | Module file      | What it contains                                                                                    | What you need in prepShinyCellPlus                                                                                                                                     | Included by data_type | Developed |
  | ----------------------- | ------------------------ | ---------------- | --------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------- | --------- |
  | `cellinfo_cellinfo`     | CellInfo vs CellInfo     | `scDRcell.R`     | 2D embedding coloured by metadata, with optional grouping or splitting by a second metadata field   | Standard ShinyCell outputs in `out_dir` (`sc1conf.rds`, `sc1meta.rds`, `sc1gene.rds`, `sc1def.rds`). Metadata columns must exist in `seurat_obj@meta.data` before prep | RNA, RNA_ATAC         | Yes       |
  | `cellinfo_geneexpr`     | CellInfo vs GeneExpr     | `scDRnum.R`      | 2D embedding with gene expression overlay (feature style plots)                                     | Standard ShinyCell outputs. `gene_mapping = TRUE` recommended if you rely on gene aliases                                                                              | RNA, RNA_ATAC         | Yes       |
  | `cellinfo3D_cellinfo3D` | CellInfo3D vs CellInfo3D | `scDRcell3D.R`   | Interactive 3D embedding coloured by metadata                                                       | `do_umap3d = TRUE` so 3D reductions exist (default from PCA). Ensure `umap3d_reductions` exists in the object                                                          | RNA, RNA_ATAC         | Yes       |
  | `cellinfo3D_geneexpr3D` | CellInfo3D vs GeneExpr3D | `scDRnum3D.R`    | Interactive 3D embedding with gene expression overlay                                               | `do_umap3d = TRUE` so 3D reductions exist                                                                                                                              | RNA, RNA_ATAC         | Yes       |
  | `genecoex`              | Gene Coexpression        | `scDRcoex.R`     | Coexpression visualisation for selected genes across cells or groups                                | Standard ShinyCell outputs. Expression must be available from prepared objects                                                                                         | RNA, RNA_ATAC         | Yes       |
  | `violin_boxplot`        | Violin / BoxPlot         | `scVioBox.R`     | Violin and boxplots for gene expression or metadata across groups                                   | Standard ShinyCell outputs                                                                                                                                             | RNA, RNA_ATAC         | Yes       |
  | `proportions`           | Cell Proportions         | `scProp.R`       | Composition summaries across groups (for example cluster proportions per sample)                    | Standard ShinyCell outputs. Requires grouping metadata such as sample and cluster                                                                                      | RNA, RNA_ATAC         | Yes       |
  | `bubble_heatmap`        | Bubble Plot / Heatmap    | `scBubbHeat.R`   | Bubble plot and heatmap summaries, typically gene sets across groups with size and colour encodings | Standard ShinyCell outputs                                                                                                                                             | RNA, RNA_ATAC         | Yes       |
  | `pseudobulk`            | Pseudobulk DE            | `scPseudobulk.R` | Pseudobulk aggregation and DE workflow based on raw counts                                          | Recommended: `do_counts_h5 = TRUE` to create `sc1counts.h5`. Requires DE packages installed (`edgeR`, `limma`)                                                         | RNA                   | Yes       |
  | `multiome_links`        | Multiome Links           | Not in list yet  | Peak to gene links and linked multiome features (RNA ATAC integration views)                        | Requires ATAC assay content and link objects prepared upstream (not generated by prepShinyCellPlus in the snippet shown)                                               | RNA_ATAC              | Not yet   |
  | `peak_browser`          | Peak Browser             | Not in list yet  | Peak level browser style views for ATAC peaks, regions, and signals                                 | Requires ATAC assay and peak level data prepared upstream                                                                                                              | RNA_ATAC              | Not yet   |
  | `eregulons_graphs`      | eRegulons Graphs         | Not in list yet  | Regulatory network and eRegulon visualisation                                                       | Requires eRegulon results and graph objects prepared upstream                                                                                                          | RNA_ATAC              | Not yet   |
  | `pseudobulk_eregulons`  | Pseudobulk eRegulons     | Not in list yet  | Pseudobulk style analysis focused on eRegulon activity                                              | Requires eRegulon matrices plus grouping metadata and counts or activity matrices                                                                                      | RNA_ATAC              | Not yet   |
  | `spatial_qc`            | Spatial QC               | Not in list yet  | Spatial quality control views (spots, metrics, filtering summaries)                                 | Requires Spatial assay and spatial metadata prepared upstream. Some parts may come from ShinyCell Spatial support                                                      | SPATIAL               | Not yet   |
  | `spatial_feature`       | Spatial Feature          | Not in list yet  | Spatial feature plots (gene expression over tissue coordinates)                                     | Requires Spatial assay and spatial coordinates prepared upstream                                                                                                       | SPATIAL               | Not yet   |
  
***

