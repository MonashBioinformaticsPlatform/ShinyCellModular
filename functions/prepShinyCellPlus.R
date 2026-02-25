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
) {
  
  .msg <- function(...) if (isTRUE(verbose)) message(...)
  .need_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
  }
  
  if (is.null(markers_file)) markers_file <- file.path(out_dir, "markergenes_lists.parquet")
  if (is.null(counts_h5_file)) counts_h5_file <- file.path(out_dir, "sc1counts.h5")
  
  ###########################################################################
  # Dependency Check
  ###########################################################################
  
  cran_pkgs <- c(
    "shiny", "shinyhelper", "shinyjs", "data.table", "Matrix", "DT",
    "magrittr", "ggplot2", "ggrepel", "hdf5r", "ggdendro", "gridExtra",
    "arrow", "rsconnect", "shinythemes", "shinydashboard", "tidyverse",
    "sortable", "plotly", "FlexDotPlot", "RColorBrewer", "ggforce"
  )
  
  bioc_pkgs <- c("limma", "edgeR")
  
  core_cran <- c("Seurat", "ShinyCell")
  cran_pkgs <- unique(c(cran_pkgs, core_cran))
  
  if (isTRUE(do_markers)) cran_pkgs <- unique(c(cran_pkgs, "presto"))
  
  missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_cran) > 0 || length(missing_bioc) > 0) {
    
    if (!isTRUE(install_missing)) {
      stop(
        paste0(
          if (length(missing_cran) > 0) paste0("Missing CRAN packages: ", paste(missing_cran, collapse = ", ")) else "",
          if (length(missing_cran) > 0 && length(missing_bioc) > 0) "\n" else "",
          if (length(missing_bioc) > 0) paste0("Missing Bioconductor packages: ", paste(missing_bioc, collapse = ", ")) else "",
          "\nSet install_missing = TRUE to install automatically."
        ),
        call. = FALSE
      )
    }
    
    if (length(missing_cran) > 0) {
      .msg("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
      install.packages(missing_cran)
    }
    
    if (length(missing_bioc) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      .msg("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
      BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
    }
    
    .msg("Dependency installation complete.")
  } else {
    .msg("All ShinyCellPlus dependencies are installed.")
  }
  
  ###########################################################################
  # Core required
  ###########################################################################
  
  .need_pkg("Seurat")
  .need_pkg("ShinyCell")
  .need_pkg("Matrix")
  
  ###########################################################################
  # Load Seurat Object
  ###########################################################################
  
  if (!is.null(seurat_rds)) {
    .msg("Loading Seurat object from: ", seurat_rds)
    seurat_obj <- readRDS(seurat_rds)
  }
  if (is.null(seurat_obj)) stop("Provide seurat_obj or seurat_rds.", call. = FALSE)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ###########################################################################
  # Ensure assay exists and set DefaultAssay
  ###########################################################################
  
  if (!assay_rna %in% Seurat::Assays(seurat_obj)) {
    stop("assay_rna not found in object: ", assay_rna, call. = FALSE)
  }
  Seurat::DefaultAssay(seurat_obj) <- assay_rna
  
  ###########################################################################
  # Optional set Idents
  ###########################################################################
  
  if (!is.null(ident_col)) {
    if (!ident_col %in% colnames(seurat_obj@meta.data)) {
      stop("ident_col not found in meta.data: ", ident_col, call. = FALSE)
    }
    Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[ident_col]]
  }
  
  ###########################################################################
  # Ensure assay keys exist
  ###########################################################################
  
  .msg("Checking assay keys")
  for (a in Seurat::Assays(seurat_obj)) {
    key <- Seurat::Key(seurat_obj[[a]])
    if (is.null(key) || !nzchar(key)) {
      Seurat::Key(seurat_obj[[a]]) <- paste0(tolower(a), "_")
      .msg("  Set Key for assay ", a, " to ", Seurat::Key(seurat_obj[[a]]))
    }
  }
  
  ###########################################################################
  # Variable features
  ###########################################################################
  
  if (isTRUE(do_variable_features)) {
    .msg("Running FindVariableFeatures on assay ", assay_rna)
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
  }
  
  ###########################################################################
  # Optional 3D UMAP
  ###########################################################################
  
  if (isTRUE(do_umap3d)) {
    .msg("Adding 3D UMAP reductions")
    for (red in umap3d_reductions) {
      if (!red %in% names(seurat_obj@reductions)) {
        .msg("  Skipping reduction ", red, " (not found)")
        next
      }
      red_name <- paste0(red, umap3d_name_suffix)
      .msg("  RunUMAP reduction=", red, " into ", red_name)
      seurat_obj <- Seurat::RunUMAP(
        seurat_obj,
        reduction = red,
        dims = umap3d_dims,
        n.components = 3,
        reduction.name = red_name
      )
    }
  }
  
  ###########################################################################
  # Optional marker genes parquet
  ###########################################################################
  
  if (isTRUE(do_markers)) {
    .need_pkg("presto")
    .need_pkg("arrow")
    
    if (file.exists(markers_file) && !isTRUE(markers_overwrite)) {
      .msg("Markers file exists, skipping (set markers_overwrite=TRUE to regenerate): ", markers_file)
    } else {
      .msg("Computing markers with presto::wilcoxauc")
      
      meta_cols <- colnames(seurat_obj@meta.data)
      resolutions <- meta_cols[grepl(markers_res_pattern, meta_cols)]
      if (!length(resolutions)) stop("No resolution columns found using pattern: ", markers_res_pattern, call. = FALSE)
      
      Seurat::DefaultAssay(seurat_obj) <- assay_rna
      expr <- Seurat::GetAssayData(seurat_obj, layer = "data")
      
      markers_list <- NULL
      for (res in resolutions) {
        .msg("  Markers for: ", res)
        clusters <- seurat_obj@meta.data[[res]]
        mk <- presto::wilcoxauc(expr, clusters)
        mk <- as.data.frame(mk)
        mk$annotation <- res
        markers_list <- if (is.null(markers_list)) mk else rbind(markers_list, mk)
      }
      
      .msg("Writing markers to: ", markers_file)
      arrow::write_parquet(markers_list, markers_file)
    }
  } else {
    .msg("Markers optional is OFF, skipping marker generation")
  }
  
  ###########################################################################
  # Optional raw counts H5 for pseudobulk (CSC)
  ###########################################################################
  
  if (isTRUE(do_counts_h5)) {
    .need_pkg("hdf5r")
    
    if (file.exists(counts_h5_file) && !isTRUE(counts_overwrite)) {
      .msg("Counts H5 exists, skipping (set counts_overwrite=TRUE to overwrite): ", counts_h5_file)
    } else {
      
      .msg("Writing sparse raw counts to H5 (CSC), file: ", counts_h5_file)
      
      counts <- Seurat::GetAssayData(seurat_obj, assay = assay_rna, layer = counts_layer)
      if (!inherits(counts, "dgCMatrix")) counts <- Matrix::as(counts, "dgCMatrix")
      
      i <- counts@i
      p <- counts@p
      x <- counts@x
      dims <- counts@Dim
      genes <- rownames(counts)
      cells <- colnames(counts)
      
      storage.mode(i) <- "integer"
      storage.mode(p) <- "integer"
      
      if (file.exists(counts_h5_file)) file.remove(counts_h5_file)
      
      h5 <- hdf5r::H5File$new(counts_h5_file, mode = "w")
      grp <- h5$create_group("counts")
      
      grp$create_dataset("i", robj = i, dtype = hdf5r::h5types$H5T_STD_I32LE, gzip_level = 4)
      grp$create_dataset("p", robj = p, dtype = hdf5r::h5types$H5T_STD_I32LE, gzip_level = 4)
      grp$create_dataset("x", robj = x, gzip_level = 4)
      grp$create_dataset("dims",  robj = as.integer(dims), dtype = hdf5r::h5types$H5T_STD_I32LE)
      grp$create_dataset("genes", robj = genes)
      grp$create_dataset("cells", robj = cells)
      
      h5$create_attr("format", "dgCMatrix_CSC_v1")
      h5$create_attr("assay", assay_rna)
      h5$create_attr("layer", counts_layer)
      
      h5$close_all()
      .msg("Counts H5 written OK")
    }
  } else {
    .msg("Counts H5 optional is OFF, skipping counts export")
  }
  
  ###########################################################################
  # ShinyCell config and app files
  ###########################################################################
  
  .msg("Creating ShinyCell config")
  scConf <- ShinyCell::createConfig(seurat_obj)
  
  if (isTRUE(do_make_app)) {
    .msg("Running makeShinyApp into: ", out_dir)
    ShinyCell::makeShinyApp(
      seurat_obj,
      scConf,
      gene.mapping = gene_mapping,
      shiny.title = shiny_title,
      shiny.dir = out_dir
    )
  } else {
    .msg("makeShinyApp optional is OFF, skipping app generation")
  }
  
  invisible(list(
    seurat_obj = seurat_obj,
    scConf = scConf,
    out_dir = out_dir,
    markers_file = if (isTRUE(do_markers)) markers_file else NULL,
    counts_h5_file = if (isTRUE(do_counts_h5)) counts_h5_file else NULL
  ))
}


read_sc1counts_h5 <- function(counts_h5_file) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("Missing package: hdf5r", call. = FALSE)
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Missing package: Matrix", call. = FALSE)
  if (!file.exists(counts_h5_file)) stop("File not found: ", counts_h5_file, call. = FALSE)
  
  h5 <- hdf5r::H5File$new(counts_h5_file, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  
  grp <- h5[["counts"]]
  i <- grp[["i"]][]
  p <- grp[["p"]][]
  x <- grp[["x"]][]
  dims <- grp[["dims"]][]
  genes <- grp[["genes"]][]
  cells <- grp[["cells"]][]
  
  m <- Matrix::sparseMatrix(
    i = as.integer(i) + 1L,
    p = as.integer(p),
    x = x,
    dims = as.integer(dims),
    repr = "C"
  )
  rownames(m) <- genes
  colnames(m) <- cells
  m
}