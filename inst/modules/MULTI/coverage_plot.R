# Coverage plot tab — genome browser style tracks, no Seurat/Signac at runtime
# id    = "coverage_plot"
# title = "Coverage Plot"
# runtime deps: Rsamtools, GenomicRanges, IRanges, RcppRoll, ggforce, patchwork, arrow, dplyr

############################################### Functions ############################################

# read Tn5 cut sites (fragment start and end) overlapping a region via tabix
# returns integer vector of positions
sc_read_cuts <- function(path, chr, start, end) {
  param  <- Rsamtools::TabixFile(path)
  region <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
  raw    <- Rsamtools::scanTabix(param, param = region)[[1]]
  if (length(raw) == 0) return(list(cuts = integer(0), cells = character(0)))
  m      <- do.call(rbind, strsplit(raw, "\t"))
  list(
    cuts  = c(as.integer(m[, 2]), as.integer(m[, 3])),  # start and end = cut sites
    cells = c(m[, 4], m[, 4])                            # cell repeated for start+end
  )
}

# compute per-group normalised coverage across a region, matching Signac's method:
#   normalised = roll_sum(raw_counts) / group_scale_factor * median_scale_factor
# where group_scale_factor = n_cells_in_group * mean(nCount_ATAC) for that group
sc_coverage <- function(frag_info, sc1meta_atac, group_col, chr, start, end, window = 100,dir_inputs) {
  
  positions <- start:end
  n_pos     <- length(positions)
  
  # compute per-group scale factors from metadata
  groups      <- unique(sc1meta_atac[[group_col]])
  groups      <- groups[!is.na(groups)]
  gsf <- vapply(groups, function(grp) {
    idx <- sc1meta_atac[[group_col]] == grp & !is.na(sc1meta_atac[[group_col]])
    n   <- sum(idx)
    if (n == 0) return(0)
    n * mean(sc1meta_atac$nCount_ATAC[idx], na.rm = TRUE)
  }, numeric(1))
  names(gsf) <- groups
  global_sf  <- median(gsf[gsf > 0])
  
  # accumulate cut counts per position per group across all fragment files
  nF <- length(frag_info)
  all_cuts <- vector("list", nF)
  frag_names <- names(frag_info)
  for (fidx in seq_len(nF)) {
    nm <- frag_names[fidx]
    # incProgress needs an active shiny::withProgress from the caller, the reactive
    # render has one, the pdf/png download handlers do not, tryCatch keeps this quiet
    # when there is no progress bar to update. this is usually the slowest part
    # (tabix scanning per fragment file), so it gets the largest share of the budget
    tryCatch(shiny::incProgress(0.6 / nF, detail = paste0("Reading fragment file ", fidx, " of ", nF)), error = function(e) NULL)
    
    fi <- frag_info[[nm]]
    # fi$path is relative to dir_inputs, resolve it there if not found as-is
    frag_path <- if (file.exists(fi$path)) fi$path else file.path(dir_inputs, fi$path)
    if (!file.exists(frag_path)){
      warning("Fragment file not found: ", fi$path, call. = FALSE)
      next
    }
    res <- tryCatch(sc_read_cuts(frag_path, chr, start, end), error = function(e) NULL)
    if (is.null(res) || length(res$cuts) == 0) next
    
    # fragment file has original barcodes; frag_info$cells maps suffixed -> original
    # invert to get original -> suffixed lookup
    orig_to_suffixed <- setNames(names(fi$cells), fi$cells)
    suffixed <- orig_to_suffixed[res$cells]          # NA if not in this sample
    keep     <- !is.na(suffixed) & suffixed %in% sc1meta_atac$cell_barcodes
    cuts_cells <- suffixed[keep]
    cuts_pos   <- res$cuts[keep]
    
    in_region  <- cuts_pos >= start & cuts_pos <= end
    cuts_cells <- cuts_cells[in_region]
    cuts_pos   <- cuts_pos[in_region]
    
    all_cuts[[fidx]] <- data.frame(pos = cuts_pos, cell = cuts_cells, stringsAsFactors = FALSE)
  }
  
  cuts_df <- do.call(rbind, Filter(Negate(is.null), all_cuts))
  if (is.null(cuts_df) || nrow(cuts_df) == 0) return(NULL)
  
  # sc1meta_atac is a data.table, so [cell_vector, col] triggers a join, not a row lookup
  # match() against cell_barcodes works the same for data.table or data.frame
  cuts_df$group <- sc1meta_atac[[group_col]][match(cuts_df$cell, sc1meta_atac$cell_barcodes)]
  cuts_df       <- cuts_df[!is.na(cuts_df$group), ]
  
  # bin and smooth per group
  nG <- length(groups)
  result_list <- vector("list", nG)
  for (gi in seq_len(nG)) {
    grp <- groups[gi]
    # incProgress needs an active shiny::withProgress from the caller, the reactive
    # render has one, the pdf/png download handlers do not, tryCatch keeps this quiet
    # when there is no progress bar to update
    tryCatch(shiny::incProgress(0.3 / nG, detail = paste0("Group ", gi, " of ", nG)), error = function(e) NULL)
    
    sub <- cuts_df[cuts_df$group == grp, ]
    if (nrow(sub) == 0) next
    
    raw <- tabulate(sub$pos - start + 1L, nbins = n_pos)
    
    # roll_sum smoothing matching Signac window
    smoothed <- RcppRoll::roll_sum(raw, n = window, fill = NA, align = "center")
    smoothed[is.na(smoothed)] <- 0
    
    sf <- if (gsf[grp] > 0) gsf[grp] else 1
    normalised <- smoothed / sf * global_sf
    
    result_list[[gi]] <- data.frame(
      pos      = positions,
      coverage = normalised,
      group    = grp,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, result_list)
}

# fast coverage using the precomputed per-chromosome parquet tables from
# prepShinyCellModular (ATAC/tiles/<chr>.parquet), only opens the one chromosome
# file needed for this region instead of scanning fragment files live, same
# normalisation as sc_coverage() so the two stay comparable
sc_coverage_tiles <- function(tiles_dir, sc1meta_atac, group_col, chr, start, end) {
  
  tile_file <- file.path(tiles_dir, paste0(chr, ".parquet"))
  if (!file.exists(tile_file)) return(NULL)
  
  # .data/.env pronouns avoid clashing "start"/"end" the columns with start/end the
  # function arguments, arrow pushes this filter down so only matching rows are read
  sub <- tryCatch(
    arrow::open_dataset(tile_file) |>
      dplyr::filter(.data$end >= .env$start, .data$start <= .env$end) |>
      dplyr::collect(),
    error = function(e) NULL
  )
  if (is.null(sub) || nrow(sub) == 0) return(NULL)
  
  sub$group <- sc1meta_atac[[group_col]][match(sub$cell, sc1meta_atac$cell_barcodes)]
  sub <- sub[!is.na(sub$group), ]
  if (nrow(sub) == 0) return(NULL)
  
  # tab stage calculation, kept here on purpose: scale factors depend on which
  # grouping/split the user picked interactively, so they cannot be precomputed at
  # prep time for every possible combination, this is cheap (metadata only, no
  # fragment or tile data involved)
  groups      <- unique(sc1meta_atac[[group_col]])
  groups      <- groups[!is.na(groups)]
  gsf <- vapply(groups, function(grp) {
    idx <- sc1meta_atac[[group_col]] == grp & !is.na(sc1meta_atac[[group_col]])
    n   <- sum(idx)
    if (n == 0) return(0)
    n * mean(sc1meta_atac$nCount_ATAC[idx], na.rm = TRUE)
  }, numeric(1))
  names(gsf) <- groups
  global_sf  <- median(gsf[gsf > 0])
  
  # bins are anchored to an arbitrary per-chromosome offset at prep time (min start
  # on that chromosome), not to the query start, so use the bin boundaries actually
  # present in the data rather than reconstructing a grid, which would not line up
  bin_grid <- unique(sub[, c("start", "end")])
  bin_grid <- bin_grid[order(bin_grid$start), ]
  
  nG <- length(groups)
  result_list <- vector("list", nG)
  for (gi in seq_len(nG)) {
    grp <- groups[gi]
    # parquet read is already fast, this is mostly for a consistent progress bar
    tryCatch(shiny::incProgress(0.3 / nG, detail = paste0("Group ", gi, " of ", nG)), error = function(e) NULL)
    
    grp_sub <- sub[sub$group == grp, ]
    
    # tab stage calculation, kept here on purpose: same reason as gsf above, this
    # depends on the interactively chosen group/split, not something prep can do
    # sum counts across cells in this group, per bin, rowsum is much faster than
    # aggregate for this, same idea as the rowsum used in prepShinyCellModular.R
    raw <- rep(0, nrow(bin_grid))
    if (nrow(grp_sub) > 0) {
      sums       <- rowsum(grp_sub$count, group = grp_sub$start, reorder = TRUE)
      bin_starts <- as.integer(rownames(sums))
      raw[match(bin_starts, bin_grid$start)] <- as.integer(sums[, 1])
    }
    
    sf <- if (gsf[grp] > 0) gsf[grp] else 1
    normalised <- raw / sf * global_sf
    
    # one row per bin, not per bp, no smoothing pass either, the bins are already
    # averaged at prep time (see atac_tile_binsize)
    result_list[[gi]] <- data.frame(
      start    = bin_grid$start,
      end      = bin_grid$end,
      coverage = normalised,
      group    = grp,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, result_list)
}

sc_plot_coverage <- function(cov_df, chr, start, end, base_size = 14) {
  if (is.null(cov_df) || nrow(cov_df) == 0)
    return(ggplot() + theme_void() +
             annotate("text", x = .5, y = .5, label = "No coverage data", colour = "grey50"))
  
  ymax   <- signif(max(cov_df$coverage, na.rm = TRUE), 2)
  cov_df$coverage <- pmin(cov_df$coverage, ymax)
  
  # preserve group order if factor
  if (!is.factor(cov_df$group))
    cov_df$group <- factor(cov_df$group, levels = unique(cov_df$group))
  
  n_grp  <- nlevels(cov_df$group)
  colors <- scales::hue_pal()(n_grp)
  names(colors) <- levels(cov_df$group)
  
  # binned data (from sc_coverage_tiles) has start/end columns, one row per bin,
  # drawn as blocks since the value really is flat within each bin. per-bp data
  # (from sc_coverage, the live tabix path) has pos, one row per base pair, drawn
  # as a continuous area like before
  if ("start" %in% names(cov_df)) {
    p <- ggplot(cov_df, aes(xmin = start, xmax = end + 1, ymin = 0, ymax = coverage, fill = group)) +
      geom_rect()
  } else {
    p <- ggplot(cov_df, aes(pos, coverage, fill = group)) +
      geom_area(stat = "identity", alpha = 1)
  }
  
  p +
    geom_hline(yintercept = 0, linewidth = 0.1) +
    facet_wrap(~ group, ncol = 1, strip.position = "left") +
    scale_fill_manual(values = colors) +
    scale_x_continuous(limits = c(start, end), expand = c(0, 0),
                       labels = scales::comma) +
    # only show 0 and the max tick, default pretty breaks (0,1000,2000...) overlap
    # once there are many facets (e.g. 13 clusters) and not enough vertical room per row
    scale_y_continuous(limits = c(0, ymax), breaks = c(0, ymax)) +
    xlab(paste0(chr, " position (bp)")) +
    ylab(paste0("Normalized signal\n(range 0 - ", ymax, ")")) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position    = "none",
      strip.background   = element_blank(),
      strip.text.y.left  = element_text(angle = 0, hjust = 1),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.line.x.bottom = element_blank(),
      panel.spacing.y    = unit(0, "lines")
    )
}

sc_plot_peaks <- function(peaks, chr, start, end, base_size = 14) {
  if (is.null(peaks)) return(NULL)
  sub <- as.data.frame(peaks[
    GenomicRanges::seqnames(peaks) == chr &
      GenomicRanges::start(peaks)   <= end  &
      GenomicRanges::end(peaks)     >= start
  ])
  if (nrow(sub) == 0) return(NULL)
  sub$start <- pmax(sub$start, start)
  sub$end   <- pmin(sub$end,   end)
  ggplot(sub) +
    geom_segment(aes(x = start, xend = end, y = 0, yend = 0),
                 linewidth = 2, colour = "dimgrey") +
    scale_x_continuous(limits = c(start, end), expand = c(0, 0)) +
    xlab(NULL) + ylab("Peaks") +
    theme_classic(base_size = base_size) +
    theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x.bottom = element_blank())
}

sc_plot_links <- function(links, chr, start, end, base_size = 14) {
  if (is.null(links) || length(links) == 0) return(NULL)
  sub <- as.data.frame(links[
    GenomicRanges::seqnames(links) == chr &
      GenomicRanges::start(links)   <= end  &
      GenomicRanges::end(links)     >= start
  ])
  if (nrow(sub) == 0) return(NULL)
  
  sub$group <- seq_len(nrow(sub))
  df <- data.frame(
    x     = c(sub$start, (sub$start + sub$end) / 2, sub$end),
    y     = c(rep(0, nrow(sub)), rep(-1, nrow(sub)), rep(0, nrow(sub))),
    group = rep(sub$group, 3),
    score = rep(sub$score, 3)
  )
  min_col <- min(0, min(df$score))
  ggplot(df) +
    ggforce::geom_bezier(aes(x = x, y = y, group = group, colour = score)) +
    geom_hline(yintercept = 0, colour = "grey") +
    scale_colour_gradient2(low = "red", mid = "grey", high = "blue",
                           limits = c(min_col, max(df$score)), n.breaks = 3) +
    scale_x_continuous(limits = c(start, end), expand = c(0, 0)) +
    xlab(NULL) + ylab("Links") +
    theme_classic(base_size = base_size) +
    theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x.bottom = element_blank())
}

sc_plot_annotation <- function(annotation, chr, start, end, base_size = 14) {
  if (is.null(annotation)) return(NULL)
  
  region_gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
  sub <- IRanges::subsetByOverlaps(annotation, region_gr)
  if (length(sub) == 0) return(NULL)
  
  sub  <- sub[sub$type == "exon"]
  if (length(sub) == 0) return(NULL)
  
  exons <- as.data.frame(sub)
  exons <- exons[!is.na(exons$gene_name), ]
  
  genes_keep <- unique(exons$gene_name)
  # get full gene bodies from the original annotation (not clipped to region)
  all_ann <- as.data.frame(
    annotation[!is.na(annotation$gene_name) &
                 annotation$gene_name %in% genes_keep &
                 annotation$type == "exon"]
  )
  
  bodies <- do.call(rbind, lapply(genes_keep, function(g) {
    rows <- all_ann[all_ann$gene_name == g, ]
    data.frame(gene_name = g, start = min(rows$start), end = max(rows$end),
               strand = rows$strand[1], stringsAsFactors = FALSE)
  }))
  
  # record overlap stacking
  bodies$dodge <- 1L
  if (nrow(bodies) > 1) {
    bodies <- bodies[order(bodies$start), ]
    row_end <- -Inf
    dodge   <- 1L
    for (i in seq_len(nrow(bodies))) {
      if (bodies$start[i] < row_end) {
        dodge <- dodge + 1L
      } else {
        dodge   <- 1L
        row_end <- bodies$end[i] + 1000L
      }
      bodies$dodge[i] <- dodge
    }
  }
  exons$dodge  <- bodies$dodge[match(exons$gene_name, bodies$gene_name)]
  bodies$position <- (pmax(bodies$start, start) + pmin(bodies$end, end)) / 2
  
  p <- ggplot() +
    geom_segment(data = exons,
                 aes(x = start, xend = end, y = dodge, yend = dodge,
                     colour = strand),
                 linewidth = 3, show.legend = FALSE) +
    geom_segment(data = bodies,
                 aes(x = pmax(start, !!start), xend = pmin(end, !!end),
                     y = dodge, yend = dodge, colour = strand),
                 linewidth = 0.5, show.legend = FALSE) +
    geom_text(data = bodies,
              aes(x = position, y = dodge + 0.3, label = gene_name),
              size = 3.5) +
    scale_colour_manual(values = c("+" = "darkblue", "-" = "darkgreen",
                                   "*" = "darkblue")) +
    scale_x_continuous(limits = c(start, end), expand = c(0, 0),
                       labels = scales::comma) +
    scale_y_continuous(limits = c(0.5, max(bodies$dodge) + 0.6)) +
    xlab(paste0(chr, " position (bp)")) + ylab("Genes") +
    theme_classic(base_size = base_size) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank())
}

sc_plot_expression <- function(gene, sc1meta, sc1gene, h5_path, group_col, base_size = 14) {
  if (is.null(gene) || !nzchar(gene) || !gene %in% names(sc1gene)) return(NULL)
  if (!file.exists(h5_path)) return(NULL)
  h5file <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
  vals   <- h5file[["grp"]][["data"]]$read(args = list(sc1gene[gene], quote(expr=)))
  vals[vals < 0] <- 0
  df <- data.frame(val = vals, group = sc1meta[[group_col]], stringsAsFactors = FALSE)
  df <- df[!is.na(df$group), ]
  # each facet should only know about its own group, otherwise facet_wrap keeps the
  # discrete y-scale shared across all facets and reserves row space for every other
  # group's level too, which is what was squeezing each violin down to a sliver
  df$group <- droplevels(factor(df$group))
  ggplot(df, aes(x = val, y = group, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
    facet_wrap(~ group, ncol = 1, strip.position = "right", scales = "free_y") +
    scale_x_continuous(position = "bottom", limits = c(0, NA)) +
    xlab(gene) + ylab(NULL) +
    theme_classic(base_size = base_size) +
    theme(legend.position    = "none",
          strip.background   = element_blank(),
          strip.text.y       = element_blank(),
          axis.text.y        = element_blank(),
          axis.ticks.y       = element_blank())
}

sc_gene_region <- function(gene, annotation) {
  if (is.null(annotation)) return(NULL)
  hits <- annotation[!is.na(annotation$gene_name) & annotation$gene_name == gene]
  if (length(hits) == 0) return(NULL)
  # some gene_name values are reused across multiple, unrelated loci (repeat families,
  # small RNA genes, pseudogenes), so restrict to the first hit's chromosome and cluster
  # nearby entries, instead of spanning min/max across every matching locus genome wide
  first_chr <- as.character(GenomicRanges::seqnames(hits)[1])
  hits      <- hits[as.character(GenomicRanges::seqnames(hits)) == first_chr]
  clusters  <- GenomicRanges::reduce(hits, min.gapwidth = 10000L)
  locus     <- IRanges::subsetByOverlaps(clusters, hits[1])
  list(chr   = first_chr,
       start = min(GenomicRanges::start(locus)),
       end   = max(GenomicRanges::end(locus)))
}

sc_parse_region <- function(s) {
  s     <- gsub(":", "-", trimws(s))
  parts <- strsplit(s, "-")[[1]]
  if (length(parts) != 3) return(NULL)
  list(chr = parts[1], start = as.integer(parts[2]), end = as.integer(parts[3]))
}

# remove x-axis from all but the last track — mirrors CombineTracks
sc_combine_tracks <- function(track_list, heights, expr_plot = NULL, widths = c(10, 2)) {
  active  <- Filter(Negate(is.null), track_list)
  h       <- heights[!sapply(track_list, is.null)]
  n       <- length(active)
  if (n == 0) return(ggplot() + theme_void())
  
  strip_x <- function(p) p + theme(
    axis.title.x       = element_blank(),
    axis.text.x        = element_blank(),
    axis.line.x.bottom = element_blank(),
    axis.ticks.x       = element_blank()
  )
  for (i in seq_len(n - 1)) active[[i]] <- strip_x(active[[i]])
  
  main <- patchwork::wrap_plots(active, ncol = 1, heights = h)
  
  if (!is.null(expr_plot)) {
    (active[[1]] + expr_plot + patchwork::plot_layout(widths = widths)) /
      patchwork::wrap_plots(active[-1], ncol = 1, heights = h[-1])
  } else {
    main
  }
}

############################################### UI ###################################################

coverage_plot_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
  tabPanel(
    HTML("Coverage Plot"),
    h4("Genome browser coverage tracks"),
    "In this tab, users can visualise Tn5 insertion frequency across genomic regions ",
    "with optional annotation, peak, link, and expression tracks.",
    br(),
    "Coverage data is precomputed and binned at 500bp resolution during data preparation, ",
    "so values are already averaged within each 500bp window rather than at single base pair precision.",
    br(), br(),
    
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        h4("Region"),
        radioButtons(ns("sc1cov_input_method"), "Input method:",
                     choices = c("Gene", "Region"), selected = "Gene", inline = TRUE),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Gene'", ns("sc1cov_input_method")),
          selectizeInput(ns("sc1cov_gene"), "Gene name:",
                         choices = NULL, selected = NULL,
                         options = list(placeholder = "Type a gene name", maxOptions = 10))
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Region'", ns("sc1cov_input_method")),
          selectInput(ns("sc1cov_chr"), "Chromosome:", choices = NULL),
          fluidRow(
            column(6, numericInput(ns("sc1cov_start"), "Start:", value = NULL)),
            column(6, numericInput(ns("sc1cov_end"),   "End:",   value = NULL))
          )
        ),
        
        numericInput(ns("sc1cov_ext_up"), "Extend upstream (bp):",   value = 1000, min = 0, step = 500),
        numericInput(ns("sc1cov_ext_dn"), "Extend downstream (bp):", value = 1000, min = 0, step = 500),
        
        hr(),
        
        h4("Grouping"),
        
        selectInput(ns("sc1cov_group"), "Group cells by:",
                    choices  = sc1conf[grp == TRUE]$UI,
                    selected = sc1def$grp1),
        
        selectInput(ns("sc1cov_splitby"), "Split tracks by (optional):",
                    choices  = c("None", sc1conf[grp == TRUE]$UI),
                    selected = "None"),
        
        hr(),
        
        h4("Tracks"),
        
        checkboxInput(ns("sc1cov_show_peaks"), "Peaks",          value = TRUE),
        checkboxInput(ns("sc1cov_show_links"), "Links",          value = TRUE),
        checkboxInput(ns("sc1cov_show_annot"), "Annotation",     value = TRUE),
        checkboxInput(ns("sc1cov_show_expr"),  "RNA expression", value = FALSE),
        
        conditionalPanel(
          condition = sprintf("input['%s']", ns("sc1cov_show_expr")),
          selectizeInput(ns("sc1cov_expr_gene"), "Gene for expression track:",
                         choices = NULL, selected = NULL,
                         options = list(placeholder = "Type a gene name", maxOptions = 10))
        ),
        
        hr(),
        actionButton(ns("sc1cov_plot_btn"), "Update plot", class = "btn btn-primary btn-lg", style = "width: 100%;"),
        # actionButton(ns("sc1cov_plot_btn"), "Update plot", class = "btn btn-primary"),
        
        radioButtons(ns("sc1cov_psz"), "Plot size:",
                     choices = c("Small", "Medium", "Large"), selected = "Medium", inline = TRUE),
        
        radioButtons(ns("sc1cov_fsz"), "Font size:",
                     choices = c("Small", "Medium", "Large"), selected = "Medium", inline = TRUE)
      ),
      
      column(
        9,
        uiOutput(ns("sc1cov_oup.ui")),
        br(),
        downloadButton(ns("sc1cov_pdf"), "Download PDF"),
        downloadButton(ns("sc1cov_png"), "Download PNG"),
        br(),
        div(style = "display:inline-block",
            numericInput(ns("sc1cov_h"), "Height (in):", width = "120px", min = 2, max = 40, value = 8, step = 0.5)),
        div(style = "display:inline-block",
            numericInput(ns("sc1cov_w"), "Width (in):",  width = "120px", min = 2, max = 40, value = 14, step = 0.5)),
        br(), br(),
        # temporary debug output for fragment path troubleshooting
        h5("Current Plot: Log info"),
        verbatimTextOutput(ns("sc1cov_debug"))
      )
    )
  )
}

############################################### Server ###############################################

coverage_plot_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs,
                                 sc1meta_atac, sc1fragmentpaths_atac, sc1annotation_atac,
                                 sc1peaks_atac, sc1links_atac, sc1atactiles = NULL) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    observe_helpers()
    
    if (!is.null(sc1meta_atac)) {
      if (!"cell_barcodes" %in% colnames(sc1meta_atac)){ sc1meta_atac$cell_barcodes <- sc1meta_atac$sampleID}
      else(message(
        "cell_barcodes column existed, we will use it as cell identifiers"
      ))
    }
    
    gene_choices <- if (!is.null(sc1annotation_atac) &&
                        "gene_name" %in% names(GenomicRanges::mcols(sc1annotation_atac))) {
      ann <- sc1annotation_atac
      # restrict to protein_coding genes when biotype is available, this keeps repeat
      # families and small RNAs (e.g. Y_RNA, 5S_rRNA) that reuse gene_name across many
      # loci out of the picker entirely, so a selected gene is always a single locus
      if ("gene_biotype" %in% names(GenomicRanges::mcols(ann)))
        ann <- ann[!is.na(ann$gene_biotype) & ann$gene_biotype == "protein_coding"]
      sort(unique(ann$gene_name[!is.na(ann$gene_name)]))
    } else names(sc1gene)
    
    # chromosome list for the Region input method, annotation first, peaks as fallback
    chr_choices <- if (!is.null(sc1annotation_atac)) {
      sort(unique(as.character(GenomicRanges::seqnames(sc1annotation_atac))))
    } else if (!is.null(sc1peaks_atac)) {
      sort(unique(as.character(GenomicRanges::seqnames(sc1peaks_atac))))
    } else character(0)
    
    # sending choices and a selected value in the same initial message often doesn't
    # render on first load with server-side selectize, onFlushed waits until the client
    # widget actually exists before we push the default gene
    session$onFlushed(function() {
      updateSelectizeInput(session, "sc1cov_gene",      choices = gene_choices,   selected = gene_choices[1], server = TRUE)
      updateSelectizeInput(session, "sc1cov_expr_gene", choices = names(sc1gene), server = TRUE)
      updateSelectInput(session, "sc1cov_chr", choices = chr_choices, selected = chr_choices[1])
    }, once = TRUE)
    
    if (!exists("pList", inherits = TRUE))
      pList <<- c(Small = "400px", Medium = "650px", Large = "900px")
    
    # per-facet row height in px, used together with n_groups_val() below so the total
    # plot height scales with how many groups (e.g. clusters) are being faceted, instead
    # of squeezing everything into the fixed pList height meant for ~6 sample groups
    if (!exists("rowPxList", inherits = TRUE))
      rowPxList <<- c(Small = 22, Medium = 32, Large = 42)
    
    # font base_size passed to theme_classic() in every track function, Medium = 14
    # keeps the old hardcoded default so nothing changes until the user touches the radio
    if (!exists("sList", inherits = TRUE))
      sList <<- c(Small = 10, Medium = 14, Large = 18)
    
    # temporary debug output for fragment path troubleshooting
    debugTxt <- reactiveVal("")
    
    # number of groups (facet rows) in the current coverage/expression plot, updated
    # inside plot_result() once cov_df is known, defaults to 1 before the first plot
    n_groups_val <- reactiveVal(1L)
    
    # eventReactive gates the coverage computation behind the Plot button, so it only
    # runs once on load and then again whenever the user clicks Plot, instead of on
    # every input change, same run-button pattern as pseudobulk.R's input$sc1e1_run
    plot_result <- eventReactive(input$sc1cov_plot_btn, {
      shiny::withProgress(message = "Loading coverage plot", value = 0, {
        
        gene <- trimws(input$sc1cov_gene %||% "")
        r    <- if (input$sc1cov_input_method %||% "Gene" == "Gene") {
          if (nzchar(gene)) sc_gene_region(gene, sc1annotation_atac) else NULL
        } else {
          chr <- input$sc1cov_chr
          st  <- input$sc1cov_start
          en  <- input$sc1cov_end
          if (!is.null(chr) && nzchar(chr) && !is.null(st) && !is.null(en) && !is.na(st) && !is.na(en))
            list(chr = chr, start = as.integer(st), end = as.integer(en))
          else NULL
        }
        # falling back to gene_choices[1] here is safe again now that this block is
        # gated behind the Plot button, it only affects the one automatic run on load
        if (is.null(r) && length(gene_choices) > 0) r <- sc_gene_region(gene_choices[1], sc1annotation_atac)
        if (is.null(r) || is.null(sc1meta_atac) || is.null(sc1fragmentpaths_atac)) {
          # temporary debug output for fragment path troubleshooting
          dbg <- paste0(
            "Early return, one of these is NULL:\n",
            "  r is.null: ", is.null(r), "\n",
            "  sc1meta_atac is.null: ", is.null(sc1meta_atac), "\n",
            "  sc1fragmentpaths_atac is.null: ", is.null(sc1fragmentpaths_atac), "\n",
            "  input method: '", input$sc1cov_input_method %||% "", "'\n",
            "  gene input: '", gene, "'\n",
            "  chr/start/end input: '", input$sc1cov_chr %||% "", "', '", input$sc1cov_start %||% "", "', '", input$sc1cov_end %||% "", "'\n",
            "  sc1annotation_atac is.null: ", is.null(sc1annotation_atac)
          )
          debugTxt(dbg)
          message(dbg)
          return(ggplot() + theme_void())
        }
        r$start <- max(1L, r$start - (input$sc1cov_ext_up %||% 1000))
        r$end   <- r$end + (input$sc1cov_ext_dn %||% 1000)
        
        grp       <- sc1conf[UI == input$sc1cov_group]$ID
        split_lbl <- input$sc1cov_splitby
        split_col <- if (!is.null(split_lbl) && split_lbl != "None") sc1conf[UI == split_lbl]$ID else NULL
        meta      <- sc1meta_atac
        meta[["__grp__"]] <- if (!is.null(split_col) && split_col %in% colnames(meta))
          paste0(meta[[split_col]], "_", meta[[grp]]) else meta[[grp]]
        
        # temporary debug output for fragment path troubleshooting
        frag_status <- vapply(names(sc1fragmentpaths_atac), function(nm) {
          fi <- sc1fragmentpaths_atac[[nm]]
          fp <- if (file.exists(fi$path)) fi$path else file.path(dir_inputs, fi$path)
          paste0("  [", nm, "] ", fp, " (exists: ", file.exists(fp), ")")
        }, character(1))
        dbg <- paste0(
          "Region: ", r$chr, ":", r$start, "-", r$end, "\n",
          "Group column: ", grp, "\n",
          #"Fragment files:\n", paste(frag_status, collapse = "\n"), "\n",
          # temporary debug output for links troubleshooting
          "sc1links_atac is.null: ", is.null(sc1links_atac), "\n",
          "sc1links_atac total length: ", if (is.null(sc1links_atac)) NA else length(sc1links_atac), "\n",
          # temporary debug output for tile matrix troubleshooting
          "sc1atactiles is.null: ", is.null(sc1atactiles), " (", if (is.null(sc1atactiles)) "using live tabix scan" else "using precomputed tile parquet", ")"
        )
        
        # temporary debug output for timing troubleshooting
        t_data0 <- Sys.time()
        cov_df <- tryCatch(
          if (!is.null(sc1atactiles))
            sc_coverage_tiles(sc1atactiles, meta, "__grp__", r$chr, r$start, r$end)
          else
            sc_coverage(sc1fragmentpaths_atac, meta, "__grp__", r$chr, r$start, r$end, 100, dir_inputs),
          error = function(e) { message("Coverage error: ", conditionMessage(e)); NULL }
        )
        t_data1 <- Sys.time()
        
        # temporary debug output for fragment path troubleshooting
        dbg <- paste0(dbg, "\n",
                      "cov_df: ", if (is.null(cov_df)) "NULL" else paste0(nrow(cov_df), " rows, groups: ", paste(unique(cov_df$group), collapse = ", ")), "\n",
                      # temporary debug output for timing troubleshooting
                      "cov_df computation time: ", round(as.numeric(difftime(t_data1, t_data0, units = "secs")), 2), "s"
        )
        debugTxt(dbg)
        message(dbg)
        
        # drives the plot height scaling in sc1cov_oup.ui, falls back to 1 if there is
        # no coverage data so the height calc below never sees a NULL/zero facet count
        n_groups_val(if (is.null(cov_df)) 1L else length(unique(cov_df$group)))
        
        shiny::incProgress(0.1, detail = "Building tracks")
        
        # temporary debug output for timing troubleshooting
        t_tracks0 <- Sys.time()
        fsz <- sList[input$sc1cov_fsz %||% "Medium"]
        track_cov   <- sc_plot_coverage(cov_df, r$chr, r$start, r$end, base_size = fsz)
        track_peaks <- if (isTRUE(input$sc1cov_show_peaks)) sc_plot_peaks(sc1peaks_atac, r$chr, r$start, r$end, base_size = fsz)           else NULL
        track_links <- if (isTRUE(input$sc1cov_show_links)) sc_plot_links(sc1links_atac, r$chr, r$start, r$end, base_size = fsz)           else NULL
        track_annot <- if (isTRUE(input$sc1cov_show_annot)) sc_plot_annotation(sc1annotation_atac, r$chr, r$start, r$end, base_size = fsz) else NULL
        track_expr  <- if (isTRUE(input$sc1cov_show_expr)) {
          g <- trimws(input$sc1cov_expr_gene %||% "")
          if (nzchar(g)) sc_plot_expression(g, sc1meta, sc1gene, file.path(dir_inputs, "RNA", "sc1gexpr.h5"), grp, base_size = fsz) else NULL
        } else NULL
        
        combined <- sc_combine_tracks(
          list(coverage = track_cov, peaks = track_peaks, links = track_links, annotation = track_annot),
          c(coverage = 10, peaks = 1, links = 2, annotation = 2),
          expr_plot = track_expr
        )
        t_tracks1 <- Sys.time()
        
        # temporary debug output for timing troubleshooting
        dbg <- paste0(dbg, "\n",
                      "track building + combine time: ", round(as.numeric(difftime(t_tracks1, t_tracks0, units = "secs")), 2), "s"
        )
        debugTxt(dbg)
        message(dbg)
        
        combined
      })
    }, ignoreNULL = FALSE)
    
    output$sc1cov_oup <- renderPlot({
      p <- plot_result()
      
      # temporary debug output for timing troubleshooting
      t_draw0 <- Sys.time()
      shiny::withProgress(message = "Rendering plot", value = 0.5, {
        print(p)
      })
      t_draw1 <- Sys.time()
      isolate({
        debugTxt(paste0(debugTxt(), "\n",
                        "draw time: ", round(as.numeric(difftime(t_draw1, t_draw0, units = "secs")), 2), "s"))
      })
    })
    
    output$sc1cov_oup.ui <- renderUI({
      # extra_px leaves room for the non-faceted tracks below the coverage/expression
      # facets (peaks, links, annotation), height_px never goes below the plain pList
      # size, it only grows past it once there are enough groups to need the room
      extra_px  <- 250
      base_px   <- as.numeric(gsub("px", "", pList[input$sc1cov_psz]))
      min_height_px <- max(base_px, n_groups_val() * rowPxList[input$sc1cov_psz] + extra_px)
      
      # sc1cov_h / sc1cov_w are the same inches used by the PDF/PNG downloads, px_per_in
      # converts them to screen pixels so the on-screen plot can be sized the same way
      # instead of only reacting to the Small/Medium/Large radio
      px_per_in <- 96
      height_px <- max(min_height_px, (input$sc1cov_h %||% 0) * px_per_in)
      width_px  <- (input$sc1cov_w %||% 0) * px_per_in
      
      plotOutput(ns("sc1cov_oup"),
                 height = paste0(height_px, "px"),
                 width  = if (width_px > 0) paste0(width_px, "px") else "100%")
    })
    
    # temporary debug output for fragment path troubleshooting
    output$sc1cov_debug <- renderText({ debugTxt() })
    
    output$sc1cov_pdf <- downloadHandler(
      filename = function() paste0("coverage_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content  = function(file) {
        gene <- trimws(input$sc1cov_gene %||% "")
        r    <- if (input$sc1cov_input_method %||% "Gene" == "Gene") {
          if (nzchar(gene)) sc_gene_region(gene, sc1annotation_atac) else NULL
        } else {
          chr <- input$sc1cov_chr
          st  <- input$sc1cov_start
          en  <- input$sc1cov_end
          if (!is.null(chr) && nzchar(chr) && !is.null(st) && !is.null(en) && !is.na(st) && !is.na(en))
            list(chr = chr, start = as.integer(st), end = as.integer(en))
          else NULL
        }
        # r should never be NULL, fall back to the first gene in the list if both gene and region inputs are empty or invalid
        if (is.null(r) && length(gene_choices) > 0) r <- sc_gene_region(gene_choices[1], sc1annotation_atac)
        if (is.null(r)) return(NULL)
        r$start <- max(1L, r$start - (input$sc1cov_ext_up %||% 1000))
        r$end   <- r$end + (input$sc1cov_ext_dn %||% 1000)
        grp       <- sc1conf[UI == input$sc1cov_group]$ID
        split_lbl <- input$sc1cov_splitby
        split_col <- if (!is.null(split_lbl) && split_lbl != "None") sc1conf[UI == split_lbl]$ID else NULL
        meta      <- sc1meta_atac
        meta[["__grp__"]] <- if (!is.null(split_col) && split_col %in% colnames(meta))
          paste0(meta[[split_col]], "_", meta[[grp]]) else meta[[grp]]
        cov_df <- tryCatch(
          if (!is.null(sc1atactiles))
            sc_coverage_tiles(sc1atactiles, meta, "__grp__", r$chr, r$start, r$end)
          else
            sc_coverage(sc1fragmentpaths_atac, meta, "__grp__", r$chr, r$start, r$end, 100, dir_inputs),
          error = function(e) NULL
        )
        fsz <- sList[input$sc1cov_fsz %||% "Medium"]
        p <- sc_combine_tracks(
          list(coverage   = sc_plot_coverage(cov_df, r$chr, r$start, r$end, base_size = fsz),
               peaks      = if (isTRUE(input$sc1cov_show_peaks)) sc_plot_peaks(sc1peaks_atac, r$chr, r$start, r$end, base_size = fsz) else NULL,
               links      = if (isTRUE(input$sc1cov_show_links)) sc_plot_links(sc1links_atac, r$chr, r$start, r$end, base_size = fsz) else NULL,
               annotation = if (isTRUE(input$sc1cov_show_annot)) sc_plot_annotation(sc1annotation_atac, r$chr, r$start, r$end, base_size = fsz) else NULL),
          c(coverage = 10, peaks = 1, links = 2, annotation = 2)
        )
        ggsave(file, plot = p, device = "pdf", height = input$sc1cov_h, width = input$sc1cov_w, useDingbats = FALSE)
      }
    )
    
    output$sc1cov_png <- downloadHandler(
      filename = function() paste0("coverage_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) {
        gene <- trimws(input$sc1cov_gene %||% "")
        r    <- if (input$sc1cov_input_method %||% "Gene" == "Gene") {
          if (nzchar(gene)) sc_gene_region(gene, sc1annotation_atac) else NULL
        } else {
          chr <- input$sc1cov_chr
          st  <- input$sc1cov_start
          en  <- input$sc1cov_end
          if (!is.null(chr) && nzchar(chr) && !is.null(st) && !is.null(en) && !is.na(st) && !is.na(en))
            list(chr = chr, start = as.integer(st), end = as.integer(en))
          else NULL
        }
        # r should never be NULL, fall back to the first gene in the list if both gene and region inputs are empty or invalid
        if (is.null(r) && length(gene_choices) > 0) r <- sc_gene_region(gene_choices[1], sc1annotation_atac)
        if (is.null(r)) return(NULL)
        r$start <- max(1L, r$start - (input$sc1cov_ext_up %||% 1000))
        r$end   <- r$end + (input$sc1cov_ext_dn %||% 1000)
        grp       <- sc1conf[UI == input$sc1cov_group]$ID
        split_lbl <- input$sc1cov_splitby
        split_col <- if (!is.null(split_lbl) && split_lbl != "None") sc1conf[UI == split_lbl]$ID else NULL
        meta      <- sc1meta_atac
        meta[["__grp__"]] <- if (!is.null(split_col) && split_col %in% colnames(meta))
          paste0(meta[[split_col]], "_", meta[[grp]]) else meta[[grp]]
        cov_df <- tryCatch(
          if (!is.null(sc1atactiles))
            sc_coverage_tiles(sc1atactiles, meta, "__grp__", r$chr, r$start, r$end)
          else
            sc_coverage(sc1fragmentpaths_atac, meta, "__grp__", r$chr, r$start, r$end, 100, dir_inputs),
          error = function(e) NULL
        )
        fsz <- sList[input$sc1cov_fsz %||% "Medium"]
        p <- sc_combine_tracks(
          list(coverage   = sc_plot_coverage(cov_df, r$chr, r$start, r$end, base_size = fsz),
               peaks      = if (isTRUE(input$sc1cov_show_peaks)) sc_plot_peaks(sc1peaks_atac, r$chr, r$start, r$end, base_size = fsz) else NULL,
               links      = if (isTRUE(input$sc1cov_show_links)) sc_plot_links(sc1links_atac, r$chr, r$start, r$end, base_size = fsz) else NULL,
               annotation = if (isTRUE(input$sc1cov_show_annot)) sc_plot_annotation(sc1annotation_atac, r$chr, r$start, r$end, base_size = fsz) else NULL),
          c(coverage = 10, peaks = 1, links = 2, annotation = 2)
        )
        ggsave(file, plot = p, device = "png", height = input$sc1cov_h, width = input$sc1cov_w)
      }
    )
  }
  )
}


############################################### Registration ##########################################

register_tab(
  id          = "coverage_plot",
  title       = "Coverage Plot (Multi)",
  ui          = coverage_plot_ui,
  server      = coverage_plot_server,
  author      = "Laura Perlaza-Jimenez",
  description = "Genome browser coverage tracks for ATAC-seq data, new tab created by MGBP",
  version     = "1.0",
  date        = "Jul 2026",
  source      = "MGBP custom",
  contact     = "laura.perlaza-jimenez@monash.edu"
)