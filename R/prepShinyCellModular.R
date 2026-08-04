#' @importFrom stats na.omit
#' @importFrom utils install.packages
NULL

#' Prepare a Seurat object for ShinyCellModular
#'
#' Takes a Seurat object from single cell experiments and prepares all files
#' needed to run a ShinyCellModular interactive Shiny app. Writes HDF5 counts,
#' ShinyCell config files, optional 3D UMAP, marker genes, and motif data to
#' \code{out_dir}.
#'
#' @param seurat_obj Seurat object. Alternative to \code{seurat_rds}.
#' @param seurat_rds Path to a \code{.rds} Seurat object. Alternative to \code{seurat_obj}.
#' @param out_dir Output directory. Default: \code{'Files_ShinyCellModular'}.
#' @param shiny_title Title for the Shiny app.
#' @param assays_selected Assay(s) to process. e.g. \code{c('RNA', 'ATAC')}. Default: \code{'RNA'}.
#' @param ident_col Column to set as Idents. Default: \code{NULL} (uses existing).
#' @param do_variable_features Run \code{FindVariableFeatures} before creating the ShinyCell config. Default: \code{TRUE}.
#' @param do_markers Compute marker genes with presto. Default: \code{FALSE}.
#' @param markers_file Path for output markers parquet file. Default: auto.
#' @param markers_overwrite Overwrite existing markers file. Default: \code{FALSE}.
#' @param markers_res_pattern Regex pattern to find resolution columns. Default: \code{'res\\\\.'}.
#' @param do_umap3d Run 3D UMAP. Default: \code{FALSE}.
#' @param umap3d_reductions Reductions to use as input for 3D UMAP. Default: \code{c('pca')}.
#' @param umap3d_dims Dims to pass to UMAP (auto-capped to available). Default: \code{1:30}.
#' @param umap3d_name_suffix Suffix appended to the 3D UMAP reduction name. Default: \code{'_umap3d'}.
#' @param do_counts_h5 Write raw counts to HDF5. Default: \code{TRUE}.
#' @param counts_h5_file Path for output H5 file. Default: auto.
#' @param counts_overwrite Overwrite existing H5 file. Default: \code{TRUE}.
#' @param counts_layer Seurat layer to use for counts. Default: \code{'counts'}.
#' @param do_make_app Run \code{makeShinyApp}. Default: \code{TRUE}.
#' @param gene_mapping Map gene names in ShinyCell. Default: \code{TRUE}.
#' 
#' @param fragments_paths Optional named list of fragment file path overrides by index.
#'   e.g. \code{list('1' = '/path/to/sample1.tsv.gz')}. If \code{NULL}, paths are copied
#'   from the original paths in the object. Default: \code{NULL}.
#' @param do_motifs Extract motifs from ATAC assay. Runs automatically when ATAC is in
#'   \code{assays_selected} and the motifs slot is populated. Set to \code{FALSE} to skip.
#'   Default: \code{'auto'}.
#' @param motifs_findmotifs Output of \code{FindMotifs()}, adds enrichment scores. Default: \code{NULL}.
#' @param motifs_overwrite Overwrite existing motif files. Default: \code{TRUE}.
#' @param atac_tile_binsize Bin width in bp for the precomputed genome-wide insertion
#'   count matrix used by the Coverage Plot tab. Smaller bins are more precise but
#'   make prep slower and the saved matrix bigger. Default: \code{500}.
#'   
#' @param custom_colors Named character vector of label -> hex color to override ShinyCell
#'   default colors. Extra names not in the data are ignored. Default: \code{NULL}.
#' @param default_genes Character vector of gene names to set as defaults in the app
#'   (e.g. \code{c('CD4', 'CD8A')}). If \code{NULL}, ShinyCell picks its own defaults.
#'   Default: \code{NULL}.
#' @param install_missing Auto-install missing packages. Default: \code{FALSE}.
#' @param verbose Print progress messages. Default: \code{TRUE}.
#'   
#' @param help Print the help message and return. Default: \code{FALSE}.
#'
#' @return Invisibly returns \code{NULL}. Writes output files to \code{out_dir}.
#'
#'#' @section Step by step:
#' \enumerate{
#'   \item \strong{Dependency check.} Builds a list of required CRAN, Bioconductor,
#'     and GitHub-only packages (\code{ShinyCell}, \code{ggseqlogo}, \code{FlexDotPlot},
#'     and \code{presto} if \code{do_markers = TRUE}; \code{Signac} is added automatically
#'     if \code{"ATAC"} is in \code{assays_selected}). Missing packages stop execution
#'     with an install hint unless \code{install_missing = TRUE}, in which case they are
#'     installed automatically.
#'   \item \strong{Load the Seurat object.} Exactly one of \code{seurat_obj} /
#'     \code{seurat_rds} must resolve to a Seurat object. If both are \code{NULL} (and
#'     \code{help} wasn't explicitly set), \code{help} is forced to \code{TRUE} and the
#'     function prints its help message and returns. Passing a file path to
#'     \code{seurat_obj} or an in-memory object to \code{seurat_rds} triggers a
#'     \code{warning()} pointing at the likely mistake.
#'   \item \strong{Validate \code{assays_selected}.} Checks every entry exists in
#'     \code{Seurat::Assays(seurat_obj)}; stops with the requested vs. available assay
#'     names otherwise.
#'   \item \strong{Ensure assay keys exist.} For any assay with a missing or empty
#'     \code{Key()}, sets one automatically (\code{tolower(assay) + "_"}) so downstream
#'     Seurat/Signac calls don't fail.
#'   \item \strong{Main loop over \code{assays_selected}.} For each \code{active_assay}
#'     (\code{DefaultAssay(seurat_obj) <- active_assay}):
#'     \itemize{
#'       \item \emph{RNA branch} (\code{active_assay == "RNA"}): creates
#'         \code{out_dir/RNA/}, optionally sets \code{Idents()} from \code{ident_col},
#'         optionally adds 3D UMAP reductions (\code{do_umap3d}, skipping any reduction
#'         not present with a message), and optionally computes markers per resolution
#'         column matched by \code{markers_res_pattern} using \code{presto::wilcoxauc}
#'         (skipping resolutions with fewer than 2 groups), writing the combined table
#'         to \code{markers_file} via \code{arrow::write_parquet()}.
#'       \item \emph{ATAC branch} (\code{active_assay == "ATAC"}): creates
#'         \code{out_dir/ATAC/}, resolves whether to run motif extraction from
#'         \code{do_motifs} (\code{"auto"} checks the \code{@@motifs} slot; otherwise
#'         \code{isTRUE(do_motifs)}), extracting PWMs and names into
#'         \code{sc1motifs.rds} and \code{sc1motifs_meta.parquet} if so, optionally
#'         joined with \code{motifs_findmotifs} enrichment scores. Always extracts
#'         static ATAC objects (\code{sc1annotation.rds}, \code{sc1peaks.rds},
#'         \code{sc1links.rds}, each skipped individually with a message if the slot is
#'         empty), and copies fragment files into
#'         \code{out_dir/ATAC/fragments/index<N>/}, using \code{fragments_paths}
#'         overrides for any file missing from the object's original path. Writes
#'         \code{sc1fragmentpaths.rds}.
#'       \item \emph{All assays}: if \code{do_make_app = TRUE}, optionally runs
#'         \code{FindVariableFeatures()}, then \code{ShinyCell::createConfig()} and
#'         \code{ShinyCell::makeShinyApp()} into \code{out_dir/<assay>/}. If
#'         \code{custom_colors} is supplied, patches matching factor-level colours into
#'         the \code{sc1conf.rds} that \code{makeShinyApp()} just wrote (since
#'         \code{makeShinyApp()} ignores any pre-patched config passed in). If
#'         \code{do_counts_h5 = TRUE}, writes the sparse counts matrix (CSC format:
#'         \code{i}, \code{p}, \code{x}, \code{dims}, \code{genes}, \code{cells}) to
#'         \code{out_dir/<assay>/sc1counts.h5} via \code{hdf5r}.
#'     }
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item \code{custom_colors} passed via \code{source("colors.R")} returns a list, not
#'     a named vector. The function detects this, extracts \code{$value} automatically,
#'     and warns you to pass \code{custom_colors = source("colors.R")$value} directly
#'     next time.
#'   \item Fragment files that can't be found anywhere (neither the original path nor a
#'     \code{fragments_paths} override) are reported with a ready-to-paste
#'     \code{fragments_paths = list(...)} snippet in the message.
#'   \item \code{install_missing = TRUE} also runs \code{setRepositories(ind = 1:3)}
#'     before \code{install.packages()} so Bioconductor dependencies (e.g. for
#'     \code{Signac}) resolve correctly.
#' }
#'
#' @examples
#' \dontrun{
#' prepShinyCellModular(
#'   seurat_rds = "seurat_object.rds",
#'   out_dir    = "my_app_files",
#'   do_umap3d  = TRUE,
#'   do_markers = TRUE
#' )
#' }
#'
#' @export
prepShinyCellModular <- function(
    seurat_obj = NULL,
    seurat_rds = NULL,
    out_dir = "ShinyCellModular_app",
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
    # motif extraction  -  runs automatically when ATAC is in assays_selected
    # and the motifs slot is populated; set to FALSE to skip
    do_motifs = "auto",
    motifs_findmotifs = NULL,   # optional: output of FindMotifs(), adds enrichment scores
    motifs_overwrite = TRUE,
    fragments_paths = NULL,     # optional named list by index to override fragment file paths
    # e.g. list("1" = "/path/to/sample1.tsv.gz", "2" = "/path/to/sample2.tsv.gz")
    # if NULL, prepShinyCellModular will try to copy from the original paths in the object
    atac_tile_binsize = 500,   # bin width in bp for the coverage plot tile matrix
    custom_colors = NULL, # add here custom colors
    default_genes=NULL,
    install_missing = FALSE,
    verbose = TRUE,
    help = FALSE
) {
  
  ###########################################################################
  # Helper Functions
  ###########################################################################

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
    if (!inherits(counts, "dgCMatrix")) counts <- methods::as(counts, "dgCMatrix")

    if (nrow(counts) == 0 || ncol(counts) == 0) {
      stop(
        "GetAssayData(assay = '", active_assay, "', layer = '", counts_layer,
        "') returned an empty matrix (", nrow(counts), " features x ", ncol(counts), " cells).\n",
        "This assay likely does not have a '", counts_layer, "' layer populated ",
        "(common for regulon/activity-score assays, which usually only carry a 'data' layer).\n",
        "Fix by passing a different counts_layer for this assay (e.g. counts_layer = 'data'), ",
        "or set do_counts_h5 = FALSE to skip H5 export for it.",
        call. = FALSE
      )
    }

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

    .msg("Creating ShinyCell config")
    scConf <- ShinyCell::createConfig(seurat_obj, maxLevels = 100)
    
    # make a folder per assay for the input files
    out_dir_path <- file.path(out_dir, active_assay)
    
    .msg("Running makeShinyApp into: ", out_dir_path)
    shinyapp_args <- list(
      seurat_obj,
      scConf,
      gex.assay    = active_assay,
      gene.mapping = gene_mapping,
      shiny.title  = shiny_title,
      shiny.dir    = out_dir_path
    )
    if (!is.null(default_genes)) shinyapp_args$default.multigene <- default_genes
    do.call(ShinyCell::makeShinyApp, shinyapp_args)
    
    # Patch custom_colors into sc1conf.rds after makeShinyApp writes it,
    # as makeShinyApp overwrites the file ignoring any pre-patched scConf
    if (!is.null(custom_colors) && length(custom_colors) > 0) {
      # custom_colors must be a named character vector; source() output is a list with $value
      if (is.list(custom_colors)) {
        warning("custom_colors is a list  -  did you use source()? Extracting $value automatically. Use custom_colors = source('colors.R')$value to avoid this warning.", call. = FALSE)
        custom_colors <- custom_colors$value
      }
      if (!is.character(custom_colors) || is.null(names(custom_colors))) {
        warning("custom_colors is not a named character vector  -  skipping color patch.", call. = FALSE)
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
  
    .resolve_per_assay <- function(x, active_assay, default) {
    if (is.list(x)) {
      if (!is.null(x[[active_assay]])) return(x[[active_assay]])
      return(default)
    }
    x
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
      .msg("No motifs slot found in assay ", active_assay, "  -  skipping motif extraction")
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
      .msg("No annotation slot found in assay ", active_assay, "  -  skipping")
    }
    
    peaks <- tryCatch(seurat_obj[[active_assay]]@ranges, error = function(e) NULL)
    if (!is.null(peaks) && length(peaks) > 0) {
      saveRDS(peaks, file.path(atac_out_dir, "sc1peaks.rds"))
      .msg("Saved peak ranges to sc1peaks.rds")
    } else {
      .msg("No ranges slot found in assay ", active_assay, "  -  skipping")
    }
    
    links <- tryCatch(Signac::Links(seurat_obj[[active_assay]]), error = function(e) NULL)
    if (!is.null(links) && length(links) > 0) {
      saveRDS(links, file.path(atac_out_dir, "sc1links.rds"))
      .msg("Saved links to sc1links.rds")
    } else {
      .msg("No links found in assay ", active_assay, "  -  skipping")
    }
    
    frags    <- tryCatch(Signac::Fragments(seurat_obj[[active_assay]]), error = function(e) NULL)
    frag_dir <- file.path(atac_out_dir, "fragments")
    
    if (is.null(frags) || length(frags) == 0) {
      .msg("No fragment files found in assay ", active_assay, "  -  skipping")
      return(invisible(NULL))
    }
    
    dir.create(frag_dir, recursive = TRUE, showWarnings = FALSE)
    orig_paths <- vapply(frags, function(f) Signac::GetFragmentData(f, slot = "path"), character(1))
    
    # check which orig_paths are missing and whether the user supplied overrides for them
    missing_orig  <- which(!file.exists(orig_paths))
    user_supplied <- vapply(seq_along(frags), function(i) !is.null(fragments_paths[[as.character(i)]]), logical(1))
    still_missing <- missing_orig[!user_supplied[missing_orig]]
    
    if (length(still_missing) > 0) {
      hint      <- paste(vapply(seq_along(frags), function(i) paste0("  \"", i, "\" = \"/path/to/fragment_", i, ".tsv.gz\""), character(1)), collapse = ",\n")
      missing_p <- paste(paste0("  [", still_missing, "] ", orig_paths[still_missing]), collapse = "\n")
      message(paste0("Fragment files not found:\n", missing_p, "\n\nRe-run with following order but replacing the paths for the location in your local computer:\nfragments_paths = list(\n", hint, "\n)"))
    }
    
    # check that user-supplied paths actually exist before proceeding
    if (!is.null(fragments_paths)) {
      bad_user_paths <- Filter(function(p) !file.exists(p), unlist(fragments_paths))
      if (length(bad_user_paths) > 0)
        stop("Supplied fragments_paths do not exist:\n",
             paste0("  ", bad_user_paths, collapse = "\n"), call. = FALSE)
    }
    
    frag_info <- lapply(seq_along(frags), function(i) {
      f         <- frags[[i]]
      src_path  <- orig_paths[i]
      src_index <- paste0(src_path, ".tbi")
      cells     <- Signac::GetFragmentData(f, slot = "cells")
      
      # each fragment gets its own index subfolder: fragments/index1/, fragments/index2/, ...
      idx_dir   <- file.path(frag_dir, paste0("index", i))
      dir.create(idx_dir, recursive = TRUE, showWarnings = FALSE)
      dst_path  <- normalizePath(file.path(idx_dir, "atac_fragments.tsv.gz"), mustWork = FALSE)
      dst_index <- paste0(dst_path, ".tbi")
      
      # path relative to dir_inputs, so containers can resolve it after the data dir is remounted
      rel_path  <- file.path("ATAC", "fragments", paste0("index", i), "atac_fragments.tsv.gz")
      
      # prefer orig_path if it exists, fall back to user-supplied
      user_path  <- fragments_paths[[as.character(i)]]
      active_src <- if (file.exists(src_path)) src_path else user_path
      active_idx <- paste0(active_src, ".tbi")
      
      if (is.null(active_src)) {
        .msg("  Fragment ", i, " not available — skipping copy")
        return(list(path = rel_path, cells = cells, index = i, copied = FALSE))
      }
      
      if (!file.exists(dst_path)) {
        file.copy(active_src, dst_path)
        .msg("  Copied fragment ", i, " to: ", dst_path)
      } else {
        .msg("  Fragment ", i, " already exists, skipping copy: ", dst_path)
      }
      if (file.exists(active_idx) && !file.exists(dst_index))
        file.copy(active_idx, dst_index)
      
      list(path = rel_path, cells = cells, index = i, copied = TRUE)
    })
    names(frag_info) <- as.character(seq_along(frags))
    
    saveRDS(frag_info, file.path(atac_out_dir, "sc1fragmentpaths.rds"))
    .msg("Saved fragment paths to sc1fragmentpaths.rds")
    invisible(frag_info)
  }
  
  # builds one binned Tn5 cut count table per chromosome (start, end, cell, count),
  # reads each fragment file once here at prep time, so coverage_plot.R only has to
  # open the one chromosome's parquet file for a region query instead of scanning
  # the fragment file fresh on every click, same idea as ArchR's per-chromosome
  # Arrow files, but stored as parquet like the rest of this project
  .extract_atac_tilematrix <- function(atac_out_dir, frag_info, bin_size = 500) {
    
    .need_pkg("data.table")
    .need_pkg("arrow")
    
    if (is.null(frag_info) || length(frag_info) == 0) {
      .msg("No fragment info available  -  skipping tile matrix")
      return(invisible(NULL))
    }
    
    
    #.msg("creating frag_dir") #debug
    frag_dir <- file.path(atac_out_dir, "fragments")
    
    # pass 1: read every fragment file once, tag each row with the suffixed cell id
    # used elsewhere in the app (same mapping sc_coverage() builds at runtime),
    # fragments that were never actually copied are skipped
    #.msg(" cuts_list") #debug
    cuts_list <- lapply(frag_info, function(fi) {
      if (!isTRUE(fi$copied)) return(NULL)
      frag_path <- file.path(frag_dir, paste0("index", fi$index), "atac_fragments.tsv.gz")
      if (!file.exists(frag_path)) return(NULL)
      
      dt <- data.table::fread(frag_path, header = FALSE,
                              select = 1:4, col.names = c("chr", "start", "end", "barcode"))
      #class(dt)
      # fi$cells: names = suffixed id used in the app, values = original barcode in the file
      orig_to_suffixed <- setNames(names(fi$cells), fi$cells)
      dt$cell <- orig_to_suffixed[dt$barcode]
      dt <- dt[!is.na(dt$cell), ]
    
      dt[, c("chr", "start", "end", "cell")]
    })
    
    
    #.msg(" cuts_all") #debug
    
    cuts_all <- data.table::rbindlist(Filter(Negate(is.null), cuts_list))
    if (nrow(cuts_all) == 0) {
      .msg("No usable fragment data found  -  skipping tile matrix")
      return(invisible(NULL))
    }
    #.msg("  tiles_dir") #debug
    tiles_dir <- file.path(atac_out_dir, "tiles")
    dir.create(tiles_dir, recursive = TRUE, showWarnings = FALSE)
    
    # one parquet file per chromosome, only cover spans that actually have fragment
    # data, so a region query at runtime only ever opens the one file it needs
    chrs <- unique(cuts_all$chr)
    .msg("  Processing ", length(chrs), "contigs/chromosomes")
    ncontigs<-0
    for (ch in chrs) {
      
      sub <- cuts_all[cuts_all$chr == ch,]
      class(cuts_all)
      lo  <- min(sub$start)
      hi  <- max(sub$end)
      bin_starts <- seq.int(lo, hi, by = bin_size)
      
      # each row contributes one Tn5 cut at start and one at end, same convention as
      # sc_read_cuts()/sc_coverage() at runtime
      cut_positions <- c(sub$start, sub$end)
      cut_cells     <- c(sub$cell, sub$cell)
      
      bin_i <- findInterval(cut_positions, bin_starts)
      keep  <- bin_i >= 1 & bin_i <= length(bin_starts)
      
      hits <- data.table::data.table(
        start = bin_starts[bin_i[keep]],
        cell  = cut_cells[keep]
      )
      
      key   <- paste(hits$start, hits$cell, sep = "\r")
      sums  <- rowsum(rep(1L, length(key)), group = key, reorder = TRUE)
      parts <- strsplit(rownames(sums), "\r", fixed = TRUE)
      counts <- data.frame(
        start = as.integer(vapply(parts, `[`, character(1), 1L)),
        cell  = vapply(parts, `[`, character(1), 2L),
        count = as.integer(sums[, 1]),
        stringsAsFactors = FALSE
      )
      counts$end <- pmin(counts$start + bin_size - 1L, hi)
      counts <- counts[, c("start", "end", "cell", "count")]
      
      arrow::write_parquet(counts, file.path(tiles_dir, paste0(ch, ".parquet")))
      ncontigs<- ncontigs +1
      .msg("  Saved tile counts for ", ch, " to tiles/", ch, ".parquet (", nrow(counts), " bin x cell entries).   Done ", ncontigs, " of ",  length(chrs))
    }
    
    .msg("Saved binned insertion counts as parquet (", length(chrs), " chromosomes, bin size ", bin_size, "bp)")
    invisible(NULL)
  }
  
  ###########################################################################
  # File path defaults
  ###########################################################################
  
 # if (is.null(markers_file))   markers_file   <- file.path(out_dir, "markergenes_lists.parquet")
#  if (is.null(counts_h5_file)) counts_h5_file <- file.path(out_dir, "sc1counts.h5")
  
  if (is.null(markers_file))   markers_file   <- file.path(out_dir, "RNA", "markergenes_lists.parquet")
  if (is.null(counts_h5_file)) counts_h5_file <- file.path(out_dir, "RNA", "sc1counts.h5")
  
  ###########################################################################
  # Dependency Check
  ###########################################################################
  
  cran_pkgs <- c(
    "shiny", "shinyhelper", "shinyjs", "data.table", "Matrix", "DT",
    "magrittr", "ggplot2", "ggrepel", "hdf5r", "ggdendro", "gridExtra",
    "arrow", "rsconnect", "shinythemes", "shinydashboard", "tidyverse",
    "sortable", "plotly", "RColorBrewer", "ggforce", "Seurat"
  )
  
 
  if ("ATAC" %in% assays_selected) cran_pkgs <- unique(c(cran_pkgs, "Signac"))
  
  
  
  # GitHub-only packages: cannot be installed via install.packages()
  for (.gh_pkg in Filter(Negate(is.null), list(
    list(pkg = "ShinyCell",  repo = "SGDDNB/ShinyCell"),
    list(pkg = "ggseqlogo",  repo = "omarwagih/ggseqlogo"),
    list(pkg = "FlexDotPlot", repo = "Simon-Leonard/FlexDotPlot"),
    if (isTRUE(do_markers))          list(pkg = "presto", repo = "immunogenomics/presto")
    
  ))) {
    if (!requireNamespace(.gh_pkg$pkg, quietly = TRUE)) {
      if (!isTRUE(install_missing)) {
        stop("Missing GitHub package: ", .gh_pkg$pkg,
             "\nInstall with: devtools::install_github('", .gh_pkg$repo, "')", call. = FALSE)
      }
      if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
      remotes::install_github(.gh_pkg$repo)
    }
  }
  # Bioconductor packages
  bioc_pkgs <- c("limma", "edgeR")
  

  
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
      setRepositories(ind = 1:3)  # needed to automatically install Bioconductor dependencies (e.g. Signac)
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
  

  if (isTRUE(help) || (is.null(seurat_obj) && is.null(seurat_rds))) {
    tools::Rd2txt(tools::Rd_db("ShinyCellModular")[["prepShinyCellModular.Rd"]])
    return(invisible(NULL))
  }
  
  
  if (!is.null(seurat_obj) && is.character(seurat_obj)) {
    warning("seurat_obj is a character path  -  did you mean to use seurat_rds instead? Passing a file path to seurat_obj will fail since it is not read from disk.", call. = FALSE)
  }
  if (!is.null(seurat_rds)) {
    if (!is.character(seurat_rds)) {
      warning("seurat_rds is not a character path  -  did you mean to use seurat_obj instead? Passing an already-loaded Seurat object to seurat_rds will fail inside readRDS().", call. = FALSE)
    } else if (!file.exists(seurat_rds)) {
      stop(
        "Cannot find seurat_rds: ", seurat_rds, "\n",
        "Current working directory: ", getwd(), "\n",
        "Check that you are running this from the correct project folder,\n",
        "or provide the full/absolute path to seurat_rds.",
        call. = FALSE
      )
    }
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

      rna_out_dir <- file.path(out_dir, "RNA")
      if (!dir.exists(rna_out_dir))
        dir.create(rna_out_dir, recursive = TRUE, showWarnings = FALSE)
      
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
          #n_dimensions  <-  ncol(Embeddings(seurat_obj, red))
          n_dimensions  <-  ncol(Embeddings(seurat_obj[[red]]))
          if(identical(n_dimensions,umap3d_dims)){umap3d_dims= umap3d_dims }else{umap3d_dims= 1:n_dimensions}
          if (length(umap3d_dims) < 3) next
          .msg("  RunUMAP reduction=", red, " into ", red_name)
          .msg(red, " with ", n_dimensions, "dimensions")
          suppressWarnings(
          seurat_obj <- Seurat::RunUMAP(
            seurat_obj,
            reduction      = red,
            dims           = umap3d_dims,
            n.components   = 3,
            reduction.name = red_name,
            umap.method    = "uwot",
            verbose        = FALSE
          ))
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
        .msg("No motifs slot found or do_motifs = FALSE  -  skipping motif extraction")
      }
      
      .msg("Extracting static ATAC objects (annotation, peaks, links, fragments)")
      #.extract_atac_static(seurat_obj, active_assay, atac_out_dir, fragments_paths)
      frag_info <- .extract_atac_static(seurat_obj, active_assay, atac_out_dir, fragments_paths)
      
      .msg("Building binned insertion matrix for fast coverage plotting")
      .extract_atac_tilematrix(atac_out_dir, frag_info, bin_size = atac_tile_binsize)
    }
    
    ###########################################################################
    # ALL ASSAYS
    ###########################################################################  
    
if (isTRUE(do_variable_features) && active_assay != "ATAC") {
  method <- if (active_assay == "RNA") "vst" else "mean.var.plot"
  .msg("Running FindVariableFeatures on assay ", active_assay, " (method: ", method, ")")
  seurat_obj <- tryCatch(
    Seurat::FindVariableFeatures(seurat_obj, selection.method = method),
    error = function(e) {
      .msg(
        "  FindVariableFeatures failed for assay ", active_assay,
        " using method '", method, "': ", conditionMessage(e),
        "\n  Skipping variable feature selection for this assay ",
        "(ShinyCell may warn about missing variable genes, but app creation will continue)."
      )
      seurat_obj
    }
  )
}
    
    if (isTRUE(do_make_app)) {
      createSCfiles(seurat_obj, active_assay, out_dir)
      cat("DEBUG: custom_colors at call time is", if(is.null(custom_colors)) "NULL" else paste(length(custom_colors), "entries"), "\n")  # debug
      
    } else {
      .msg("makeShinyApp optional is OFF, skipping app generation")
    }
    
    if (isTRUE(do_counts_h5)) {
      .need_pkg("hdf5r")
      h5_path <- 
        file.path(out_dir, active_assay, "sc1counts.h5")
        this_layer <- .resolve_per_assay(counts_layer, active_assay, "data")
      createcountsh5(seurat_obj, h5_path, counts_overwrite, this_layer, active_assay)
      #createcountsh5(seurat_obj, h5_path, counts_overwrite, counts_layer, active_assay)
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
