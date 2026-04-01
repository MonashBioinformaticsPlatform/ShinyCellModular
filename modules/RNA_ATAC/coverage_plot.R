# Coverage plot tab — genome browser style tracks, no Seurat/Signac at runtime
# id    = "coverage_plot"
# title = "Coverage Plot"
# runtime deps: Rsamtools, GenomicRanges, IRanges, RcppRoll, ggforce, patchwork

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
sc_coverage <- function(frag_info, sc1meta_atac, group_col, chr, start, end, window = 100) {
  
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
  all_cuts <- lapply(names(frag_info), function(nm) {
    fi <- frag_info[[nm]]
    if (!file.exists(fi$path)) {
      warning("Fragment file not found: ", fi$path, call. = FALSE)
      return(NULL)
    }
    res <- tryCatch(sc_read_cuts(fi$path, chr, start, end), error = function(e) NULL)
    if (is.null(res) || length(res$cuts) == 0) return(NULL)
    
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
    
    data.frame(pos = cuts_pos, cell = cuts_cells, stringsAsFactors = FALSE)
  })
  
  cuts_df <- do.call(rbind, Filter(Negate(is.null), all_cuts))
  if (is.null(cuts_df) || nrow(cuts_df) == 0) return(NULL)
  
  cuts_df$group <- sc1meta_atac[cuts_df$cell, group_col]
  cuts_df       <- cuts_df[!is.na(cuts_df$group), ]
  
  # bin and smooth per group
  do.call(rbind, lapply(groups, function(grp) {
    sub <- cuts_df[cuts_df$group == grp, ]
    if (nrow(sub) == 0) return(NULL)
    
    raw <- tabulate(sub$pos - start + 1L, nbins = n_pos)
    
    # roll_sum smoothing matching Signac window
    smoothed <- RcppRoll::roll_sum(raw, n = window, fill = NA, align = "center")
    smoothed[is.na(smoothed)] <- 0
    
    sf <- if (gsf[grp] > 0) gsf[grp] else 1
    normalised <- smoothed / sf * global_sf
    
    data.frame(
      pos      = positions,
      coverage = normalised,
      group    = grp,
      stringsAsFactors = FALSE
    )
  }))
}

sc_plot_coverage <- function(cov_df, chr, start, end) {
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
  
  ggplot(cov_df, aes(pos, coverage, fill = group)) +
    geom_area(stat = "identity", alpha = 1) +
    geom_hline(yintercept = 0, linewidth = 0.1) +
    facet_wrap(~ group, ncol = 1, strip.position = "left") +
    scale_fill_manual(values = colors) +
    scale_x_continuous(limits = c(start, end), expand = c(0, 0),
                       labels = scales::comma) +
    ylim(c(0, ymax)) +
    xlab(paste0(chr, " position (bp)")) +
    ylab(paste0("Normalized signal\n(range 0 - ", ymax, ")")) +
    theme_classic() +
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

sc_plot_peaks <- function(peaks, chr, start, end) {
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
    theme_classic() +
    theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x.bottom = element_blank())
}

sc_plot_links <- function(links, chr, start, end) {
  if (is.null(links) || length(links) == 0) return(NULL)
  sub <- as.data.frame(links[
    GenomicRanges::seqnames(links) == chr &
      GenomicRanges::start(links)   <= end  &
      GenomicRanges::end(links)     >= start
  ])
  if (nrow(sub) == 0) return(NULL)
  sub <- sub[sub$start >= start & sub$end <= end, ]
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
    theme_classic() +
    theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x.bottom = element_blank())
}

sc_plot_annotation <- function(annotation, chr, start, end) {
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
              size = 2.5) +
    scale_colour_manual(values = c("+" = "darkblue", "-" = "darkgreen",
                                   "*" = "darkblue")) +
    scale_x_continuous(limits = c(start, end), expand = c(0, 0),
                       labels = scales::comma) +
    scale_y_continuous(limits = c(0.5, max(bodies$dodge) + 0.6)) +
    xlab(paste0(chr, " position (bp)")) + ylab("Genes") +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank())
}

sc_plot_expression <- function(gene, sc1meta, sc1gene, h5_path, group_col) {
  if (is.null(gene) || !nzchar(gene) || !gene %in% names(sc1gene)) return(NULL)
  if (!file.exists(h5_path)) return(NULL)
  h5file <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
  vals   <- h5file[["grp"]][["data"]]$read(args = list(sc1gene[gene], quote(expr=)))
  vals[vals < 0] <- 0
  df <- data.frame(val = vals, group = sc1meta[[group_col]], stringsAsFactors = FALSE)
  df <- df[!is.na(df$group), ]
  ggplot(df, aes(x = val, y = group, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
    facet_wrap(~ group, ncol = 1, strip.position = "right") +
    scale_x_continuous(position = "bottom", limits = c(0, NA)) +
    xlab(gene) + ylab(NULL) +
    theme_classic() +
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
  list(chr   = as.character(GenomicRanges::seqnames(hits)[1]),
       start = min(GenomicRanges::start(hits)),
       end   = max(GenomicRanges::end(hits)))
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
    br(), br(),
    
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        
        h4("Region"),
        
        selectizeInput(ns("sc1cov_gene"), "Gene name:",
                       choices = NULL, selected = NULL,
                       options = list(placeholder = "Type a gene name", maxOptions = 10)),
        
        tags$p(style = "text-align:center; color:#888; font-size:12px;", "— or —"),
        
        textInput(ns("sc1cov_region"), "Manual region (chr-start-end):",
                  placeholder = "e.g. chr1-791818-1020120"),
        
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
        
        numericInput(ns("sc1cov_window"), "Smoothing window (bp):", value = 100, min = 10, step = 50),
        
        radioButtons(ns("sc1cov_psz"), "Plot size:",
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
            numericInput(ns("sc1cov_h"), "Height (in):", width = "120px", min = 2, max = 40, value = 12, step = 0.5)),
        div(style = "display:inline-block",
            numericInput(ns("sc1cov_w"), "Width (in):",  width = "120px", min = 2, max = 40, value = 10, step = 0.5))
      )
    )
  )
}

############################################### Server ###############################################

coverage_plot_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs,
                                 sc1meta_atac, sc1fragmentpaths, sc1annotation,
                                 sc1peaks, sc1links) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    observe_helpers()
    
    if (!is.null(sc1meta_atac)) {
      if (!"cell_barcodes" %in% colnames(sc1meta_atac)){ sc1meta_atac$cell_barcodes <- sc1meta_atac$sampleID}
      else(message(
        "cell_barcodes column existed, we will use it as cell identifiers"
      ))
    }
    
    gene_choices <- if (!is.null(sc1annotation) &&
                        "gene_name" %in% names(GenomicRanges::mcols(sc1annotation))) {
      sort(unique(sc1annotation$gene_name[!is.na(sc1annotation$gene_name)]))
    } else names(sc1gene)
    
    updateSelectizeInput(session, "sc1cov_gene",      choices = gene_choices,   server = TRUE)
    updateSelectizeInput(session, "sc1cov_expr_gene", choices = names(sc1gene), server = TRUE)
    
    if (!exists("pList", inherits = TRUE))
      pList <<- c(Small = "400px", Medium = "650px", Large = "900px")
    
    output$sc1cov_oup <- renderPlot({
      gene <- trimws(input$sc1cov_gene %||% "")
      r    <- if (nzchar(gene)) {sc_gene_region(gene, sc1annotation)} else{sc_parse_region(trimws(input$sc1cov_region %||% ""))}
      if (is.null(r) || is.null(sc1meta_atac) || is.null(sc1fragmentpaths))
        return(ggplot() + theme_void())
      r$start <- max(1L, r$start - (input$sc1cov_ext_up %||% 1000))
      r$end   <- r$end + (input$sc1cov_ext_dn %||% 1000)
      
      grp       <- sc1conf[UI == input$sc1cov_group]$ID
      split_lbl <- input$sc1cov_splitby
      split_col <- if (!is.null(split_lbl) && split_lbl != "None") sc1conf[UI == split_lbl]$ID else NULL
      meta      <- sc1meta_atac
      meta[["__grp__"]] <- if (!is.null(split_col) && split_col %in% colnames(meta))
        paste0(meta[[split_col]], "_", meta[[grp]]) else meta[[grp]]
      
      cov_df <- tryCatch(
        sc_coverage(sc1fragmentpaths, meta, "__grp__", r$chr, r$start, r$end, input$sc1cov_window %||% 100),
        error = function(e) { message("Coverage error: ", conditionMessage(e)); NULL }
      )
      
      track_cov   <- sc_plot_coverage(cov_df, r$chr, r$start, r$end)
      track_peaks <- if (isTRUE(input$sc1cov_show_peaks)) sc_plot_peaks(sc1peaks, r$chr, r$start, r$end)           else NULL
      track_links <- if (isTRUE(input$sc1cov_show_links)) sc_plot_links(sc1links, r$chr, r$start, r$end)           else NULL
      track_annot <- if (isTRUE(input$sc1cov_show_annot)) sc_plot_annotation(sc1annotation, r$chr, r$start, r$end) else NULL
      track_expr  <- if (isTRUE(input$sc1cov_show_expr)) {
        g <- trimws(input$sc1cov_expr_gene %||% "")
        if (nzchar(g)) sc_plot_expression(g, sc1meta, sc1gene, file.path(dir_inputs, "sc1gexpr.h5"), grp) else NULL
      } else NULL
      
      sc_combine_tracks(
        list(coverage = track_cov, peaks = track_peaks, links = track_links, annotation = track_annot),
        c(coverage = 10, peaks = 1, links = 2, annotation = 2),
        expr_plot = track_expr
      )
    })
    
    output$sc1cov_oup.ui <- renderUI({
      plotOutput(ns("sc1cov_oup"), height = pList[input$sc1cov_psz])
    })
    
    output$sc1cov_pdf <- downloadHandler(
      filename = function() paste0("coverage_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content  = function(file) {
        gene <- trimws(input$sc1cov_gene %||% "")
        r    <- if (nzchar(gene)) sc_gene_region(gene, sc1annotation)
        else sc_parse_region(trimws(input$sc1cov_region %||% ""))
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
          sc_coverage(sc1fragmentpaths, meta, "__grp__", r$chr, r$start, r$end, input$sc1cov_window %||% 100),
          error = function(e) NULL
        )
        p <- sc_combine_tracks(
          list(coverage   = sc_plot_coverage(cov_df, r$chr, r$start, r$end),
               peaks      = if (isTRUE(input$sc1cov_show_peaks)) sc_plot_peaks(sc1peaks, r$chr, r$start, r$end) else NULL,
               links      = if (isTRUE(input$sc1cov_show_links)) sc_plot_links(sc1links, r$chr, r$start, r$end) else NULL,
               annotation = if (isTRUE(input$sc1cov_show_annot)) sc_plot_annotation(sc1annotation, r$chr, r$start, r$end) else NULL),
          c(coverage = 10, peaks = 1, links = 2, annotation = 2)
        )
        ggsave(file, plot = p, device = "pdf", height = input$sc1cov_h, width = input$sc1cov_w, useDingbats = FALSE)
      }
    )
    
    output$sc1cov_png <- downloadHandler(
      filename = function() paste0("coverage_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) {
        gene <- trimws(input$sc1cov_gene %||% "")
        r    <- if (nzchar(gene)) sc_gene_region(gene, sc1annotation)
        else sc_parse_region(trimws(input$sc1cov_region %||% ""))
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
          sc_coverage(sc1fragmentpaths, meta, "__grp__", r$chr, r$start, r$end, input$sc1cov_window %||% 100),
          error = function(e) NULL
        )
        p <- sc_combine_tracks(
          list(coverage   = sc_plot_coverage(cov_df, r$chr, r$start, r$end),
               peaks      = if (isTRUE(input$sc1cov_show_peaks)) sc_plot_peaks(sc1peaks, r$chr, r$start, r$end) else NULL,
               links      = if (isTRUE(input$sc1cov_show_links)) sc_plot_links(sc1links, r$chr, r$start, r$end) else NULL,
               annotation = if (isTRUE(input$sc1cov_show_annot)) sc_plot_annotation(sc1annotation, r$chr, r$start, r$end) else NULL),
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
      id     = "coverage_plot",
      title  = "Coverage Plot",
      ui     = coverage_plot_ui,
      server = coverage_plot_server
    )