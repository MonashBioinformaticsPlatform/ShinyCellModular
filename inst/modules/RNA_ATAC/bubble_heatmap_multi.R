############################################### Functions ###########################################

# Get gene list
scGeneList <- function(inp, inpGene) {
  
  if (is.null(inp) || is.na(inp) || !nzchar(inp)) {
    return(data.table::data.table(gene = character(0), present = logical(0)))
  }
  
  toks <- unlist(strsplit(inp, "[,;\\n\\r]+", perl = TRUE), use.names = FALSE)
  toks <- trimws(toks)
  toks <- toks[nzchar(toks)]
  
  geneList <- data.table::data.table(
    gene = unique(toks),
    present = TRUE
  )
  
  geneList[!gene %in% names(inpGene), present := FALSE]
  geneList[]
}

# Build a dendrogram grob from ggdendro data (orientation: "row" = right side, "col" = top)
# Row dendrogram sits to the RIGHT of the main plot; gene labels stay on the main plot y-axis.
# Tree grows left-to-right (root on the left, leaves aligning with the plot rows).
scDendroGrob <- function(hcData, orientation = c("row", "col"), base_size = 11) {
  orientation <- match.arg(orientation)
  seg  <- hcData$segments
  labs <- hcData$labels
  if (orientation == "row") {
    gg <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = seg,
                            ggplot2::aes(x = y, y = x, xend = yend, yend = xend)) +
      ggplot2::scale_y_continuous(breaks = labs$x, labels = labs$label,
                                  expand = c(0, 0.5)) +
      ggplot2::scale_x_continuous(expand = c(0.05, 0)) +   # root left, leaves right
      ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 4, 0, 0))  # gap on right only
  } else {
    gg <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = seg,
                            ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
      ggplot2::scale_x_continuous(breaks = labs$x, labels = labs$label,
                                  expand = c(0.05, 0)) +
      ggplot2::scale_y_continuous(expand = c(0.05, 0)) +
      ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(plot.margin = ggplot2::margin(4, 0, 0, 0))
  }
  ggplot2::ggplotGrob(gg)
}

# Plot gene expression bubbleplot / heatmap
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt,
                       inpsub1, inpsub2, inpH5, inpGene,
                       inpScl, inpRow, inpCol,
                       inpcols, inpfsz,
                       inpFacet = "None",
                       inpEngine = c("Classic ggplot", "FlexDotPlot"),
                       x_order = NULL,
                       # ── Dendrograms ──
                       inpDendRow = FALSE,
                       inpDendCol = FALSE,
                       save = FALSE) {
  
  inpEngine <- match.arg(inpEngine)
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  geneList <- scGeneList(inp, inpGene)
  geneList <- geneList[present == TRUE]
  
  #shiny::validate(shiny::need(nrow(geneList) <= 50, "More than 50 genes to plot. Please reduce the gene list."))
  shiny::validate(shiny::need(nrow(geneList) > 1, "Please input at least 2 genes to plot."))
  
  # Determine whether faceting is active
  useFacet <- !is.null(inpFacet) && nzchar(inpFacet) && inpFacet != "None"
  
  h5file <- H5File$new(inpH5, mode = "r")
  on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
  h5data <- h5file[["grp"]][["data"]]
  
  
  ggData <- data.table::data.table()
  for (iGene in geneList$gene) {
    tmp <- inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) <- c("sampleID", "sub")
    tmp$grpBy <- inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName <- iGene
    tmp$val <- h5data$read(args = list(inpGene[iGene], quote(expr=)))
    
    # Add facet column if active
    if (useFacet) {
      tmp$facetBy <- inpMeta[[inpConf[UI == inpFacet]$ID]]
    }
    
    ggData <- data.table::rbindlist(list(ggData, tmp))
  }
  
  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2]
  }
  
  shiny::validate(shiny::need(data.table::uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot."))
  
  # Build the grouping columns for aggregation
  grpCols <- c("geneName", "grpBy")
  if (useFacet) grpCols <- c(grpCols, "facetBy")
  
  ggData$val <- expm1(ggData$val)
  ggData <- ggData[, .(
    val  = mean(val),
    prop = sum(val > 0) / length(sampleID)
  ), by = grpCols]
  ggData$val <- log1p(ggData$val)
  
  colRange <- range(ggData$val)
  if (isTRUE(inpScl)) {
    ggData[, val := as.numeric(scale(val)), keyby = "geneName"]
    colRange <- c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))
  }
  
  # Clustering: use the full (non-faceted) matrix for gene/column ordering
  ggMatData <- ggData[, .(val = mean(val)), by = c("geneName", "grpBy")]
  ggMat <- data.table::dcast.data.table(ggMatData, geneName ~ grpBy, value.var = "val")
  tmp <- ggMat$geneName
  ggMat <- as.matrix(ggMat[, -1])
  rownames(ggMat) <- tmp
  
  if (isTRUE(inpRow)) {
    hcRow <- ggdendro::dendro_data(as.dendrogram(stats::hclust(stats::dist(ggMat))))
    ggData$geneName <- factor(ggData$geneName, levels = hcRow$labels$label)
  } else {
    hcRow <- NULL   # needed for dendrogram guard below
    ggData$geneName <- factor(ggData$geneName, levels = rev(geneList$gene))
  }
  
  if (isTRUE(inpCol)) {
    hcCol <- ggdendro::dendro_data(as.dendrogram(stats::hclust(stats::dist(t(ggMat)))))
    ggData$grpBy <- factor(ggData$grpBy, levels = hcCol$labels$label)
  } else {
    hcCol <- NULL   # needed for dendrogram guard below
  }
  
  # Apply manual X ordering (only if not clustering columns)
  if (!isTRUE(inpCol) && !is.null(x_order) && length(x_order) > 0) {
    x_order <- x_order[x_order %in% unique(as.character(ggData$grpBy))]
    if (length(x_order) > 0) {
      ggData$grpBy <- factor(ggData$grpBy, levels = x_order)
    }
  }
  
  if (inpPlt == "Bubbleplot") {
    
    if (inpEngine == "FlexDotPlot") {
      
      shiny::validate(
        shiny::need(requireNamespace("FlexDotPlot", quietly = TRUE),
                    "FlexDotPlot engine selected but the FlexDotPlot package is not installed.")
      )
      
      # FlexDotPlot does not natively support faceting;
      # fall back to Classic ggplot when a facet variable is active
      if (useFacet) {
        shiny::showNotification(
          "FlexDotPlot does not support faceting. Falling back to Classic ggplot.",
          type = "warning", duration = 5
        )
      } else {
        dot_df <- as.data.frame(ggData[, .(grpBy, geneName, val, prop)])
        dot_df$grpBy <- as.factor(dot_df$grpBy)
        dot_df$geneName <- as.factor(dot_df$geneName)
        
        res <- FlexDotPlot::dot_plot(
          data.to.plot = dot_df,
          size_var     = "prop",
          col_var      = "val",
          scale.by     = "radius",
          cols.use     = cList[[inpcols]],
          plot.legend  = TRUE,
          do.return    = TRUE,
          do.plot      = FALSE,
          x.lab.pos    = "bottom",
          y.lab.pos    = "left"
        )
        
        if ("dot.plot" %in% names(res)) return(res$dot.plot)
        if ("plot" %in% names(res)) return(res$plot)
        
        gg_idx <- which(vapply(res, inherits, logical(1), what = "ggplot"))
        if (length(gg_idx) > 0) return(res[[gg_idx[1]]])
        
        stop("FlexDotPlot returned an unexpected object structure.")
      }
    }
    
    ggOut <- ggplot2::ggplot(ggData, ggplot2::aes(grpBy, geneName, color = val, size = prop)) +
      ggplot2::geom_point() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0.5)) +
      ggplot2::scale_size_continuous(
        "proportion",
        range  = c(0, 8),
        limits = c(0, 1),
        breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)
      ) +
      ggplot2::scale_color_gradientn(
        "expression",
        limits  = colRange,
        colours = cList[[inpcols]]
      ) +
      ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = 15)) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), legend.box = "vertical")
    
  } else {
    
    ggOut <- ggplot2::ggplot(ggData, ggplot2::aes(grpBy, geneName, fill = val)) +
      ggplot2::geom_tile() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0.5)) +
      ggplot2::scale_fill_gradientn(
        "expression",
        limits  = colRange,
        colours = cList[[inpcols]]
      ) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 15)) +
      ggplot2::theme(axis.title = ggplot2::element_blank())
  }
  
  # Apply faceting if a facet variable is selected
  if (useFacet) {
    ggOut <- ggOut +
      ggplot2::facet_grid(. ~ facetBy, scales = "free_x", space = "free_x") +
      ggplot2::theme(
        strip.text       = ggplot2::element_text(size = sList[inpfsz] * 0.9, face = "bold"),
        strip.background = ggplot2::element_blank()
      )
  }
  
  ggLeg <- g_legend(ggOut)
  ggOut <- ggOut + ggplot2::theme(legend.position = "none")
  
  # ── Assemble layout: optionally attach dendrograms around the main plot ──
  # Strategy: use gtable surgery so dendrograms span ONLY the panel rows/cols of the
  # main plot gtable, not the axis-label rows. This prevents the dendrogram from being
  # allocated extra height/width that belongs to tick labels or axis titles.
  showDendRow <- isTRUE(inpDendRow) && !is.null(hcRow)
  showDendCol <- isTRUE(inpDendCol) && !is.null(hcCol)
  
  if (!showDendRow && !showDendCol) {
    # No dendrograms: original two-grob layout unchanged
    ggOut <- gridExtra::arrangeGrob(ggOut, ggLeg,
                                    heights = c(7, 2),
                                    layout_matrix = rbind(c(1), c(2)))
  } else {
    # Convert main plot to gtable so we can inspect row/col structure
    gt <- ggplot2::ggplotGrob(ggOut)
    
    # Find the panel row(s) and col(s) in the gtable (named "panel" in the layout)
    panel_rows <- unique(gt$layout$t[gt$layout$name == "panel"])
    panel_cols <- unique(gt$layout$l[gt$layout$name == "panel"])
    
    if (showDendRow) {
      dendR <- scDendroGrob(hcRow, "row", base_size = sList[inpfsz])
      # Add a new column to the right, width proportional to the plot null unit
      gt <- gtable::gtable_add_cols(gt, widths = grid::unit(0.2, "null"), pos = -1)
      # Insert dendrogram grob spanning ONLY the panel rows
      gt <- gtable::gtable_add_grob(gt, dendR,
                                    t = min(panel_rows), b = max(panel_rows),
                                    l = ncol(gt), r = ncol(gt),
                                    clip = "off", name = "dendro-row")
    }
    
    if (showDendCol) {
      dendC <- scDendroGrob(hcCol, "col", base_size = sList[inpfsz])
      # Add a new row above the gtable
      gt <- gtable::gtable_add_rows(gt, heights = grid::unit(0.15, "null"), pos = 0)
      # Panel cols may have shifted by the column addition above; re-query
      panel_cols_now <- unique(gt$layout$l[gt$layout$name == "panel"])
      # Insert dendrogram grob spanning ONLY the panel cols, in the new top row
      gt <- gtable::gtable_add_grob(gt, dendC,
                                    t = 1, b = 1,
                                    l = min(panel_cols_now), r = max(panel_cols_now),
                                    clip = "off", name = "dendro-col")
    }
    
    # Re-attach legend below using arrangeGrob
    ggOut <- gridExtra::arrangeGrob(gt, ggLeg,
                                    heights = c(7, 2),
                                    layout_matrix = rbind(c(1), c(2)))
  }
  
  ggOut
}

################################################# UI #################################################

scBubbHeat_ui <- function(id, sc1conf, sc1def) {
  
  ns <- shiny::NS(id)
  
  tabPanel(
    HTML("Bubbleplot / Heatmap"),
    h4("Gene expression bubbleplot / heatmap"),
    "In this tab, users can visualise the gene expression patterns of ",
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(),
    "The normalised expression are averaged, log-transformed and then plotted.",
    br(), br(),
    
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        
        textAreaInput(
          ns("sc1d1inp"),
          HTML("List of gene names <br />
               ( separated <br />
               by , or ; or newline):"),
          height = "200px",
          value = paste0(sc1def$genes, collapse = ", ")
        ),
        
        selectInput(
          ns("sc1d1grp"), "Group by:",
          choices = sc1conf[grp == TRUE]$UI,
          selected = sc1conf[grp == TRUE]$UI[1]
        ),
        
        # ── NEW: Facet variable selector ──
        selectInput(
          ns("sc1d1facet"), "Facet by (split panels):",
          choices = c("None", sc1conf[grp == TRUE]$UI),
          selected = "None"
        ),
        
        radioButtons(
          ns("sc1d1plt"), "Plot type:",
          choices = c("Bubbleplot", "Heatmap"),
          selected = "Bubbleplot", inline = TRUE
        ),
        
        checkboxInput(ns("sc1d1scl"), "Scale gene expression", value = TRUE),
        checkboxInput(ns("sc1d1row"), "Cluster rows (genes)", value = TRUE),
        checkboxInput(ns("sc1d1col"), "Cluster columns (samples)", value = FALSE),
        
        # ── NEW: Dendrogram toggles (visible only when clustering is active) ──
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("sc1d1row")),
          checkboxInput(ns("sc1d1dendRow"), "Show row dendrogram", value = FALSE)
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("sc1d1col")),
          checkboxInput(ns("sc1d1dendCol"), "Show column dendrogram", value = FALSE)
        ),
        
        # ── Order X axis (drag-and-drop) ──
        actionButton(ns("sc1d1togOrderX"), "Order X axis"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1d1togOrderX")),
          h5("Drag to reorder X axis groups"),
          uiOutput(ns("sc1d1xorder.ui")),
          actionButton(ns("sc1d1xorder_reset"), "Reset to default", class = "btn btn-primary")
        ),
        
        br(),
        actionButton(ns("sc1d1togL"), "Filter Cells"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1d1togL")),
          selectInput(
            ns("sc1d1sub1"), "Cell information to subset:",
            choices = sc1conf[grp == TRUE]$UI,
            selected = sc1def$grp1
          ),
          uiOutput(ns("sc1d1sub1.ui")),
          actionButton(ns("sc1d1sub1all"), "Select all groups", class = "btn btn-primary"),
          actionButton(ns("sc1d1sub1non"), "Deselect all groups", class = "btn btn-primary")
        ),
        
        br(), br(),
        actionButton(ns("sc1d1tog"), "Customize plot"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1d1tog")),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Bubbleplot'", ns("sc1d1plt")),
            radioButtons(
              ns("sc1d1engine"),
              "Bubbleplot engine:",
              choices  = c("Classic ggplot", "FlexDotPlot"),
              selected = "Classic ggplot",
              inline   = TRUE
            )
          ),
          
          radioButtons(
            ns("sc1d1cols"), "Colour scheme:",
            choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
            selected = "Blue-Yellow-Red"
          ),
          radioButtons(
            ns("sc1d1psz"), "Plot size:",
            choices = c("Small", "Medium", "Large"),
            selected = "Medium", inline = TRUE
          ),
          
          # ── NEW: Manual plot size override (useful for large gene lists) ──
          checkboxInput(ns("sc1d1manualSize"), "Set custom plot size", value = FALSE),
          conditionalPanel(
            condition = sprintf("input['%s'] == true", ns("sc1d1manualSize")),
            numericInput(ns("sc1d1manualH"), "Plot height (px):",
                         value = 650, min = 200, max = 4000, step = 50, width = "160px"),
            numericInput(ns("sc1d1manualW"), "Plot width (px):",
                         value = 900, min = 200, max = 4000, step = 50, width = "160px")
          ),
          
          radioButtons(
            ns("sc1d1fsz"), "Font size:",
            choices = c("Small", "Medium", "Large"),
            selected = "Medium", inline = TRUE
          )
        )
      ),
      
      column(
        9,
        h4(htmlOutput(ns("sc1d1oupTxt"))),
        uiOutput(ns("sc1d1oup.ui")),
        downloadButton(ns("sc1d1oup.pdf"), "Download PDF"),
        downloadButton(ns("sc1d1oup.png"), "Download PNG"),
        br(),
        div(
          style = "display:inline-block",
          numericInput(ns("sc1d1oup.h"), "PDF / PNG height:", width = "138px",
                       min = 4, max = 20, value = 10, step = 0.5)
        ),
        div(
          style = "display:inline-block",
          numericInput(ns("sc1d1oup.w"), "PDF / PNG width:", width = "138px",
                       min = 4, max = 20, value = 10, step = 0.5)
        )
      )
    )
  )
}

############################################## Server ################################################

scBubbHeat_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  shiny::moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    observe_helpers()
    optCrt <- "{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
    
    updateSelectizeInput(session, "sc1a1inp2", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene1, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1a3inp1", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene1, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1a3inp2", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene2, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1b2inp1", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene1, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1b2inp2", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene2, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1c1inp2", server = TRUE,
                         choices = c(sc1conf[is.na(fID)]$UI, names(sc1gene)),
                         selected = sc1conf[is.na(fID)]$UI[1],
                         options = list(
                           maxOptions = length(sc1conf[is.na(fID)]$UI) + 3,
                           create = TRUE, persist = TRUE, render = I(optCrt)
                         ))
    
    if (!exists("pList3", inherits = TRUE)) {
      pList3 <<- c(Small = "450px", Medium = "650px", Large = "850px")
    }
    
    # ── X axis ordering ──
    x_levels_default <- reactive({
      req(input$sc1d1grp)
      x_id <- sc1conf[UI == input$sc1d1grp]$ID
      req(length(x_id) == 1, !is.na(x_id), x_id != "")
      levs <- levels(sc1meta[[x_id]])
      if (is.null(levs)) levs <- sort(unique(as.character(sc1meta[[x_id]])))
      levs
    })
    
    x_order <- reactiveVal(NULL)
    
    observeEvent(x_levels_default(), {
      x_order(x_levels_default())
    }, ignoreInit = FALSE)
    
    output$sc1d1xorder.ui <- renderUI({
      req(x_levels_default())
      ord <- x_order()
      if (is.null(ord) || length(ord) == 0) ord <- x_levels_default()
      ord <- ord[ord %in% x_levels_default()]
      if (length(ord) == 0) ord <- x_levels_default()
      sortable::rank_list(
        text = "X axis group order",
        labels = ord,
        input_id = ns("sc1d1xorder_rank")
      )
    })
    
    observeEvent(input$sc1d1xorder_rank, {
      req(input$sc1d1xorder_rank)
      x_order(input$sc1d1xorder_rank)
    })
    
    observeEvent(input$sc1d1xorder_reset, {
      x_order(x_levels_default())
    })
    
    x_order_final <- reactive({
      levs <- x_levels_default()
      ord <- x_order()
      if (is.null(ord) || length(ord) == 0) return(levs)
      ord <- ord[ord %in% levs]
      if (length(ord) == 0) levs else ord
    })
    
    # ── Filter cells UI ──
    output$sc1d1sub1.ui <- renderUI({
      req(input$sc1d1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(ns("sc1d1sub2"), "Select which cells to show",
                         inline = TRUE, choices = sub, selected = sub)
    })
    
    observeEvent(input$sc1d1sub1non, {
      req(input$sc1d1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1d1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = NULL, inline = TRUE)
    })
    
    observeEvent(input$sc1d1sub1all, {
      req(input$sc1d1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1d1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = sub, inline = TRUE)
    })
    
    # ── Output text ──
    output$sc1d1oupTxt <- renderUI({
      geneList <- scGeneList(input$sc1d1inp, sc1gene)
      
      #if (nrow(geneList) > 50) {
      # return(HTML("More than 50 input genes. Please reduce the gene list."))
      #}
      
      ok <- geneList[present == TRUE]
      notok <- geneList[present == FALSE]
      
      oup <- paste0(nrow(ok), " genes OK and will be plotted")
      if (nrow(notok) > 0) {
        oup <- paste0(
          oup, "<br/>",
          nrow(notok), " genes not found (",
          paste0(notok$gene, collapse = ", "),
          ")"
        )
      }
      HTML(oup)
    })
    
    # ── Render plot ──
    output$sc1d1oup <- renderPlot({
      
      engine <- if (is.null(input$sc1d1engine) || !nzchar(input$sc1d1engine)) {
        "Classic ggplot"
      } else {
        input$sc1d1engine
      }
      
      facet_var <- if (is.null(input$sc1d1facet)) "None" else input$sc1d1facet
      
      p <- scBubbHeat(
        sc1conf, sc1meta,
        input$sc1d1inp, input$sc1d1grp, input$sc1d1plt,
        input$sc1d1sub1, input$sc1d1sub2,
        file.path(dir_inputs, "sc1gexpr.h5"),
        sc1gene,
        input$sc1d1scl, input$sc1d1row, input$sc1d1col,
        input$sc1d1cols, input$sc1d1fsz,
        inpFacet = facet_var,
        inpEngine = engine,
        x_order = x_order_final(),
        # ── Dendrograms ──
        inpDendRow = isTRUE(input$sc1d1dendRow),
        inpDendCol = isTRUE(input$sc1d1dendCol)
      )
      
      if (inherits(p, "grob") || inherits(p, "gtable")) {
        grid::grid.newpage()
        grid::grid.draw(p)
      } else {
        print(p)
      }
    })
    
    # ── Plot panel size: preset or manual override ──
    output$sc1d1oup.ui <- renderUI({
      if (isTRUE(input$sc1d1manualSize) &&
          !is.null(input$sc1d1manualH) && !is.null(input$sc1d1manualW)) {
        plotOutput(ns("sc1d1oup"),
                   height = paste0(input$sc1d1manualH, "px"),
                   width  = paste0(input$sc1d1manualW, "px"))
      } else {
        plotOutput(ns("sc1d1oup"), height = pList3[input$sc1d1psz])
      }
    })
    
    # ── Download PDF ──
    output$sc1d1oup.pdf <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1d1plt, "_", input$sc1d1grp, ".pdf")
      },
      content = function(file) {
        
        engine <- if (is.null(input$sc1d1engine) || !nzchar(input$sc1d1engine)) {
          "Classic ggplot"
        } else {
          input$sc1d1engine
        }
        
        facet_var <- if (is.null(input$sc1d1facet)) "None" else input$sc1d1facet
        
        ggplot2::ggsave(
          file, device = "pdf",
          height = input$sc1d1oup.h, width = input$sc1d1oup.w,
          plot = scBubbHeat(
            sc1conf, sc1meta,
            input$sc1d1inp, input$sc1d1grp, input$sc1d1plt,
            input$sc1d1sub1, input$sc1d1sub2,
            file.path(dir_inputs, "sc1gexpr.h5"),
            sc1gene,
            input$sc1d1scl, input$sc1d1row, input$sc1d1col,
            input$sc1d1cols, input$sc1d1fsz,
            inpFacet = facet_var,
            inpEngine = engine,
            x_order = x_order_final(),
            # ── Dendrograms ──
            inpDendRow = isTRUE(input$sc1d1dendRow),
            inpDendCol = isTRUE(input$sc1d1dendCol),
            save = TRUE
          )
        )
      }
    )
    
    # ── Download PNG ──
    output$sc1d1oup.png <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1d1plt, "_", input$sc1d1grp, ".png")
      },
      content = function(file) {
        
        engine <- if (is.null(input$sc1d1engine) || !nzchar(input$sc1d1engine)) {
          "Classic ggplot"
        } else {
          input$sc1d1engine
        }
        
        facet_var <- if (is.null(input$sc1d1facet)) "None" else input$sc1d1facet
        
        ggplot2::ggsave(
          file, device = "png",
          height = input$sc1d1oup.h, width = input$sc1d1oup.w,
          plot = scBubbHeat(
            sc1conf, sc1meta,
            input$sc1d1inp, input$sc1d1grp, input$sc1d1plt,
            input$sc1d1sub1, input$sc1d1sub2,
            file.path(dir_inputs, "sc1gexpr.h5"),
            sc1gene,
            input$sc1d1scl, input$sc1d1row, input$sc1d1col,
            input$sc1d1cols, input$sc1d1fsz,
            inpFacet = facet_var,
            inpEngine = engine,
            x_order = x_order_final(),
            # ── Dendrograms ──
            inpDendRow = isTRUE(input$sc1d1dendRow),
            inpDendCol = isTRUE(input$sc1d1dendCol),
            save = TRUE
          )
        )
      }
    )
    
  })
}

############################################### Registration #################################################

register_tab(
  id     = "bubble_heatmap_multi",
  title  = "Bubble Plot / Heatmap ",
  ui     = scBubbHeat_ui,
  server = scBubbHeat_server
)