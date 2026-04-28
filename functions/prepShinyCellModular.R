prepShinyCellModular <- function(
    seurat_obj = NULL,
    seurat_rds = NULL,
    out_dir = "Files_ShinyCell",
    shiny_title = "ShinyCellModular Intermediate",
    assays_selected = "RNA", #c("RNA","ATAC","regulon","chromvar") 
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
    # motif extraction — runs automatically when ATAC is in assays_selected
    # and the motifs slot is populated; set to FALSE to skip
    do_motifs = "auto",
    motifs_findmotifs = NULL,   # optional: output of FindMotifs(), adds enrichment scores
    motifs_overwrite = TRUE,
    fragments_paths = NULL,     # optional named list by index to override fragment file paths
    # e.g. list("1" = "/path/to/sample1.tsv.gz", "2" = "/path/to/sample2.tsv.gz")
    # if NULL, prepShinyCellModular will try to copy from the original paths in the object
    custom_colors = NULL, # add here custom colors
    help = FALSE
) {
  
  ###########################################################################
  # Helper Functions
  ###########################################################################
  helpMessage <- paste0(
    "prepShinyCellModular() — Prepare a Seurat object for ShinyCellModular\n",
    "\n",
    "ARGUMENTS\n",
    "  seurat_obj            Seurat object (alternative to seurat_rds)\n",
    "  seurat_rds            Path to .rds Seurat object (alternative to seurat_obj)\n",
    "  out_dir               Output directory. Default: 'Files_ShinyCell'\n",
    "  shiny_title           Title for the Shiny app\n",
    "  assays_selected       Assay(s) to process. e.g. c('RNA','ATAC')\n",
    "  ident_col             Column to set as Idents. Default: NULL (uses existing)\n",
    "  do_variable_features  Run FindVariableFeatures before createConfig. Default: TRUE\n",
    "  do_markers            Compute marker genes with presto. Default: FALSE\n",
    "  markers_file          Path for output markers parquet. Default: auto\n",
    "  markers_overwrite     Overwrite existing markers file. Default: FALSE\n",
    "  markers_res_pattern   Regex pattern to find resolution columns. Default: 'res\\\\.'\n",
    "  do_umap3d             Run 3D UMAP. Default: FALSE\n",
    "  umap3d_reductions     Reductions to use as input for 3D UMAP. Default: c('pca')\n",
    "  umap3d_dims           Dims to use (auto-capped to available). Default: 1:30\n",
    "  umap3d_name_suffix    Suffix for 3D UMAP reduction name. Default: '_umap3d'\n",
    "  do_counts_h5          Write raw counts to HDF5. Default: TRUE\n",
    "  counts_h5_file        Path for output H5 file. Default: auto\n",
    "  counts_overwrite      Overwrite existing H5 file. Default: TRUE\n",
    "  counts_layer          Seurat layer to use for counts. Default: 'counts'\n",
    "  do_make_app           Run makeShinyApp. Default: TRUE\n",
    "  gene_mapping          Map gene names in ShinyCell. Default: TRUE\n",
    "  custom_colors         Named character vector of label -> hex color to override\n",
    "                        ShinyCell default colors. Extra names not in data are ignored\n",
    "  install_missing       Auto-install missing packages. Default: FALSE\n",
    "  verbose               Print progress messages. Default: TRUE\n",
    "  do_motifs             Extract motifs from ATAC assay. Default: 'auto'\n",
    "  motifs_findmotifs     Output of FindMotifs() to add enrichment scores. Default: NULL\n",
    "  motifs_overwrite      Overwrite existing motif files. Default: TRUE\n",
    "  fragments_paths       Named list of fragment file path overrides by index.\n",
    "                        e.g. list('1' = '/path/to/sample1.tsv.gz'). Default: NULL\n",
    "  help                  Print this help message. Default: FALSE\n"
  )
  
  .msg <- function(...) if (isTRUE(verbose)) message(...)
  .need_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
  }
  
  createcountsh5 <- function(seurat_obj, counts_h5_file, counts_overwrite, counts_layer, active_assay) {
    if (file.exists(counts_h5_file) && !isTRUE(counts_overwrite)) {
      .msg("Counts H5 exists, skipping (set counts_overwrite=TRUE to overwrite): ", counts_h5_file)
    } else {
      .msg("Writing sparse raw counts to H5 (CSC), file: ", counts_h5_file)
      counts <- Seurat::GetAssayData(seurat_obj, assay = active_assay, layer = counts_layer)
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
      grp$create_dataset("i",     robj = i,                dtype = hdf5r::h5types$H5T_STD_I32LE, gzip_level = 4)
      grp$create_dataset("p",     robj = p,                dtype = hdf5r::h5types$H5T_STD_I32LE, gzip_level = 4)
      grp$create_dataset("x",     robj = x,                gzip_level = 4)
      grp$create_dataset("dims",  robj = as.integer(dims), dtype = hdf5r::h5types$H5T_STD_I32LE)
      grp$create_dataset("genes", robj = genes)
      grp$create_dataset("cells", robj = cells)
      h5$create_attr("format", "dgCMatrix_CSC_v1")
      h5$create_attr("assay",  active_assay)
      h5$create_attr("layer",  counts_layer)
      h5$close_all()
      .msg("Counts H5 written OK")
    }
  }
  
  createSCfiles <- function(seurat_obj, active_assay, out_dir) {
    if (isTRUE(do_variable_features)) {
      .msg("Running FindVariableFeatures on assay ", active_assay)
      seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
    }
    .msg("Creating ShinyCell config")
    scConf <- ShinyCell::createConfig(seurat_obj)
    
    if (active_assay == "RNA") {
      out_dir_path <- out_dir
    } else {
      out_dir_path <- file.path(out_dir, active_assay)
    }
    .msg("Running makeShinyApp into: ", out_dir_path)
    ShinyCell::makeShinyApp(
      seurat_obj,
      scConf,
      gex.assay    = active_assay,
      gene.mapping = gene_mapping,
      shiny.title  = shiny_title,
      shiny.dir    = out_dir_path
    )
    
    # Patch custom_colors into sc1conf.rds after makeShinyApp writes it,
    # as makeShinyApp overwrites the file ignoring any pre-patched scConf
    if (!is.null(custom_colors) && length(custom_colors) > 0) {
      # custom_colors must be a named character vector; source() output is a list with $value
      if (is.list(custom_colors)) {
        warning("custom_colors is a list — did you use source()? Extracting $value automatically. Use custom_colors = source('colors.R')$value to avoid this warning.", call. = FALSE)
        custom_colors <- custom_colors$value
      }
      if (!is.character(custom_colors) || is.null(names(custom_colors))) {
        warning("custom_colors is not a named character vector — skipping color patch.", call. = FALSE)
      } else {
        conf_path <- file.path(out_dir_path, "sc1conf.rds")
        .msg("Patching custom_colors into: ", conf_path)
        scConf <- readRDS(conf_path)
        for (i in seq_len(nrow(scConf))) {
          if (is.na(scConf$fCL[i])) next
          lvls        <- strsplit(scConf$fID[i], "\\|")[[1]]
          cols        <- strsplit(scConf$fCL[i], "\\|")[[1]]
          if (length(lvls) != length(cols)) next
          names(cols) <- lvls
          matched     <- intersect(lvls, names(custom_colors))
          if (!length(matched)) next
          cols[matched]  <- custom_colors[matched]
          scConf$fCL[i]  <- paste(unname(cols), collapse = "|")
          .msg("  Patched colors for column: ", scConf$UI[i],
               " (", length(matched), "/", length(lvls), " levels matched)")
        }
        saveRDS(scConf, conf_path)
        .msg("Patched sc1conf.rds saved to: ", conf_path)
      }
    }
  }
  
  .extract_motifs <- function(seurat_obj, active_assay, atac_out_dir,
                              findmotifs_df, overwrite) {
    
    motifs_rds  <- file.path(atac_out_dir, "sc1motifs.rds")
    motifs_parq <- file.path(atac_out_dir, "sc1motifs_meta.parquet")
    
    if (file.exists(motifs_rds) && file.exists(motifs_parq) && !isTRUE(overwrite)) {
      .msg("Motif files exist, skipping (set motifs_overwrite=TRUE to overwrite)")
      return(invisible(NULL))
    }
    
    motif_obj <- tryCatch(
      seurat_obj[[active_assay]]@motifs,
      error = function(e) NULL
    )
    
    if (is.null(motif_obj)) {
      .msg("No motifs slot found in assay ", active_assay, " — skipping motif extraction")
      return(invisible(NULL))
    }
    
    .need_pkg("arrow")
    
    pwm_list <- motif_obj@pwm
    .msg("Extracting ", length(pwm_list), " PWMs")
    
    saveRDS(pwm_list, motifs_rds)
    .msg("Saved PWM list to: ", motifs_rds)
    
    motif_ids   <- names(pwm_list)
    motif_names <- motif_ids
    tf_names    <- unlist(motif_obj@motif.names[motif_ids])
    
    meta <- data.frame(
      motif_id   = motif_ids,
      motif_name = motif_names,
      tf_name    = tf_names,
      stringsAsFactors = FALSE
    )
    
    if (!is.null(findmotifs_df)) {
      .msg("Joining FindMotifs enrichment scores")
      fm <- as.data.frame(findmotifs_df)
      fm$motif_id <- rownames(fm)
      
      enr_col  <- intersect(c("fold.enrichment", "enrichment"), names(fm))[1]
      pval_col <- intersect(c("p.value", "pvalue", "pval"),      names(fm))[1]
      padj_col <- intersect(c("p.adjust", "padj", "adj.pval"),   names(fm))[1]
      
      fm_sub <- fm[, na.omit(c("motif_id", enr_col, pval_col, padj_col)), drop = FALSE]
      names(fm_sub)[names(fm_sub) == enr_col]  <- "enrichment_score"
      names(fm_sub)[names(fm_sub) == pval_col] <- "pval"
      names(fm_sub)[names(fm_sub) == padj_col] <- "padj"
      
      meta <- merge(meta, fm_sub, by = "motif_id", all.x = TRUE)
    } else {
      meta$enrichment_score <- NA_real_
      meta$pval             <- NA_real_
      meta$padj             <- NA_real_
    }
    
    arrow::write_parquet(meta, motifs_parq)
    .msg("Saved motif metadata to: ", motifs_parq)
    invisible(list(pwm_list = pwm_list, meta = meta))
  }
  
  .extract_atac_static <- function(seurat_obj, active_assay, atac_out_dir, fragments_paths = NULL) {
    
    .need_pkg("Signac")
    
    annotation <- tryCatch(seurat_obj[[active_assay]]@annotation, error = function(e) NULL)
    if (!is.null(annotation) && length(annotation) > 0) {
      saveRDS(annotation, file.path(atac_out_dir, "sc1annotation.rds"))
      .msg("Saved annotation to sc1annotation.rds")
    } else {
      .msg("No annotation slot found in assay ", active_assay, " — skipping")
    }
    
    peaks <- tryCatch(seurat_obj[[active_assay]]@ranges, error = function(e) NULL)
    if (!is.null(peaks) && length(peaks) > 0) {
      saveRDS(peaks, file.path(atac_out_dir, "sc1peaks.rds"))
      .msg("Saved peak ranges to sc1peaks.rds")
    } else {
      .msg("No ranges slot found in assay ", active_assay, " — skipping")
    }
    
    links <- tryCatch(Signac::Links(seurat_obj[[active_assay]]), error = function(e) NULL)
    if (!is.null(links) && length(links) > 0) {
      saveRDS(links, file.path(atac_out_dir, "sc1links.rds"))
      .msg("Saved links to sc1links.rds")
    } else {
      .msg("No links found in assay ", active_assay, " — skipping")
    }
    
    frags    <- tryCatch(Signac::Fragments(seurat_obj[[active_assay]]), error = function(e) NULL)
    frag_dir <- file.path(atac_out_dir, "fragments")
    
    if (is.null(frags) || length(frags) == 0) {
      .msg("No fragment files found in assay ", active_assay, " — skipping")
      return(invisible(NULL))
    }
    
    dir.create(frag_dir, recursive = TRUE, showWarnings = FALSE)
    orig_paths <- vapply(frags, function(f) Signac::GetFragmentData(f, slot = "path"), character(1))
    
    orig_paths  <- vapply(frags, function(f) Signac::GetFragmentData(f, slot = "path"), character(1))
    user_supplied <- vapply(seq_along(frags), function(i) !is.null(fragments_paths[[as.character(i)]]), logical(1))
    missing_idx <- which(!file.exists(orig_paths) & !user_supplied)
    if (length(missing_idx) > 0) {
      hint      <- paste(vapply(seq_along(frags), function(i) paste0("  \"", i, "\" = \"", orig_paths[i], "\""), character(1)), collapse = ",
")
      missing_p <- paste(paste0("  [", missing_idx, "] ", orig_paths[missing_idx]), collapse = "
")
      message(paste0("Fragment files not found:
", missing_p, "

Re-run with:
fragments_paths = list(
", hint, "
)"))
    }
    
    frag_info <- lapply(seq_along(frags), function(i) {
      f         <- frags[[i]]
      src_path  <- orig_paths[i]
      src_index <- paste0(src_path, ".tbi")
      cells     <- Signac::GetFragmentData(f, slot = "cells")
      
      user_path <- fragments_paths[[as.character(i)]]
      if (!is.null(user_path)) {
        .msg("  Using supplied path for fragment ", i, ": ", user_path)
        return(list(path = user_path, cells = cells, copied = FALSE))
      }
      
      fname     <- paste0("fragment_", i, "_", basename(src_path))
      dst_path  <- file.path(frag_dir, fname)
      dst_index <- paste0(dst_path, ".tbi")
      
      if (file.exists(src_path)) {
        if (!file.exists(dst_path)) {
          file.copy(src_path, dst_path)
          .msg("  Copied fragment file: ", fname)
        } else {
          .msg("  Fragment file already exists, skipping copy: ", fname)
        }
        if (file.exists(src_index) && !file.exists(dst_index))
          file.copy(src_index, dst_index)
        list(path = dst_path, cells = cells, copied = TRUE)
      } else {
        list(path = dst_path, cells = cells, copied = FALSE)
      }
    })
    names(frag_info) <- as.character(seq_along(frags))
    
    saveRDS(frag_info, file.path(atac_out_dir, "sc1fragmentpaths.rds"))
    .msg("Saved fragment paths to sc1fragmentpaths.rds")
    invisible(frag_info)
  }
  
  ###########################################################################
  # File path defaults
  ###########################################################################
  
  if (is.null(markers_file))   markers_file   <- file.path(out_dir, "markergenes_lists.parquet")
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
    .msg("All ShinyCellModular dependencies are installed.")
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
  
  if (is.null(seurat_obj) && is.null(seurat_rds) && !isTRUE(help)) {
    help <- TRUE
  }
  if (isTRUE(help)) {
    message(helpMessage)
    return(invisible(NULL))
  }
  
  if (!is.null(seurat_rds)) {
    .msg("Loading Seurat object from: ", seurat_rds)
    seurat_obj <- readRDS(seurat_rds)
  }
  if (is.null(seurat_obj)) stop("Provide seurat_obj or seurat_rds. To get more information about the function prep", call. = FALSE)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ###########################################################################
  # Validate assays_selected
  ###########################################################################
  
  if (!all(assays_selected %in% Seurat::Assays(seurat_obj))) {
    assay_missing <- assays_selected[!(assays_selected %in% Seurat::Assays(seurat_obj))]
    print(paste("assays_selected:", assays_selected))
    print(paste("assays found in the seurat object:", Seurat::Assays(seurat_obj)))
    stop("assays_selected not found in object: ", assay_missing, call. = FALSE)
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
  # Main loop over assays
  ###########################################################################
  
  for (active_assay in assays_selected) {
    
    Seurat::DefaultAssay(seurat_obj) <- active_assay
    
    
    ###########################################################################
    # RNA ASSAY
    ###########################################################################  
    if (active_assay == "RNA") {
      
      if (!is.null(ident_col)) {
        if (!ident_col %in% colnames(seurat_obj@meta.data))
          stop("ident_col not found in meta.data: ", ident_col, call. = FALSE)
        Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[ident_col]]
      }
      
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
            reduction      = red,
            dims           = umap3d_dims,
            n.components   = 3,
            reduction.name = red_name
          )
        }
      }
      
      if (isTRUE(do_markers)) {
        .need_pkg("presto")
        .need_pkg("arrow")
        if (file.exists(markers_file) && !isTRUE(markers_overwrite)) {
          .msg("Markers file exists, skipping (set markers_overwrite=TRUE to regenerate): ", markers_file)
        } else {
          .msg("Computing markers with presto::wilcoxauc")
          meta_cols   <- colnames(seurat_obj@meta.data)
          resolutions <- meta_cols[grepl(markers_res_pattern, meta_cols)]
          if (!length(resolutions))
            stop("No resolution columns found using pattern: ", markers_res_pattern, call. = FALSE)
          Seurat::DefaultAssay(seurat_obj) <- active_assay
          expr <- Seurat::GetAssayData(seurat_obj, layer = "data")
          markers_list <- NULL
          
            for (res in resolutions) {
              clusters <- seurat_obj@meta.data[[res]]
              if (length(unique(clusters)) < 2) {
                .msg("  Skipping markers for: ", res, " (fewer than 2 unique groups)")
                next
              }
              .msg("  Markers for: ", res)
              clusters <- seurat_obj@meta.data[[res]]
              mk       <- presto::wilcoxauc(expr, clusters)
              mk       <- as.data.frame(mk)
              mk$annotation <- res
              markers_list <- if (is.null(markers_list)) mk else rbind(markers_list, mk)
            }
          .msg("Writing markers to: ", markers_file)
          arrow::write_parquet(markers_list, markers_file)
         
        }
      } else {
        .msg("Markers optional is OFF, skipping marker generation")
      }
    }
    
    
    ###########################################################################
    # ATAC ASSAY
    ###########################################################################  
    if (active_assay == "ATAC") {
      
      atac_out_dir <- file.path(out_dir, "ATAC")
      if (!dir.exists(atac_out_dir))
        dir.create(atac_out_dir, recursive = TRUE, showWarnings = FALSE)
      
      run_motifs <- if (identical(do_motifs, "auto")) {
        !is.null(tryCatch(seurat_obj[[active_assay]]@motifs, error = function(e) NULL))
      } else {
        isTRUE(do_motifs)
      }
      
      if (run_motifs) {
        .msg("Extracting motifs from assay: ", active_assay)
        .extract_motifs(
          seurat_obj    = seurat_obj,
          active_assay  = active_assay,
          atac_out_dir  = atac_out_dir,
          findmotifs_df = motifs_findmotifs,
          overwrite     = motifs_overwrite
        )
      } else {
        .msg("No motifs slot found or do_motifs = FALSE — skipping motif extraction")
      }
      
      .msg("Extracting static ATAC objects (annotation, peaks, links, fragments)")
      .extract_atac_static(seurat_obj, active_assay, atac_out_dir, fragments_paths)
    }
    
    ###########################################################################
    # ALL ASSAYS
    ###########################################################################  
    
    
    if (isTRUE(do_make_app)) {
      createSCfiles(seurat_obj, active_assay, out_dir)
      cat("DEBUG: custom_colors at call time is", if(is.null(custom_colors)) "NULL" else paste(length(custom_colors), "entries"), "\n")  # debug
      
    } else {
      .msg("makeShinyApp optional is OFF, skipping app generation")
    }
    
    if (isTRUE(do_counts_h5)) {
      .need_pkg("hdf5r")
      h5_path <- if (active_assay == "RNA") {
        counts_h5_file
      } else {
        file.path(out_dir, active_assay, "sc1counts.h5")
      }
      createcountsh5(seurat_obj, h5_path, counts_overwrite, counts_layer, active_assay)
    } else {
      .msg("Counts H5 optional is OFF, skipping counts export")
    }
  }
  
  ###########################################################################
  # Return
  ###########################################################################
  
  invisible(list(
    seurat_obj     = seurat_obj,
    out_dir        = out_dir,
    markers_file   = if (isTRUE(do_markers)) markers_file else NULL,
    counts_h5_file = if (isTRUE(do_counts_h5)) counts_h5_file else NULL
  ))

  #######################################################################
  #
  # to print help and details on the function:
  #            prepShinyCellModular(help=TRUE) 
  #
  ########################################################################

  
}