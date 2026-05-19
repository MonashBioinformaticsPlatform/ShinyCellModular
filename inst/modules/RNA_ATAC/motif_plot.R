# Motif logo plot tab
# id    = "motif_plot"
# title = "Motif Plot"
# requires: sc1motifs.rds and sc1motifs_meta.parquet in dir_inputs/ATAC/

############################################### Functions ############################################

sc_motif_plot <- function(pwm_list, meta_df, motif_ids, ncol, show_enrichment) {
  
  if (length(motif_ids) == 0) return(ggplot() + theme_void())
  
  .need_ggseqlogo <- function() {
    if (!requireNamespace("ggseqlogo", quietly = TRUE))
      stop("Package ggseqlogo is required for motif plots. Install with: install.packages('ggseqlogo')", call. = FALSE)
  }
  .need_ggseqlogo()
  pwms <- pwm_list[motif_ids]
  labels <- setNames(
    paste0(meta_df$motif_name[match(motif_ids, meta_df$motif_id)],
           "\n", motif_ids),
    motif_ids
  )
  
  logo_plots <- lapply(motif_ids, function(mid) {
    ggseqlogo::ggseqlogo(pwms[[mid]], method = "prob") +
      ggtitle(labels[mid]) +
      theme(
        plot.title   = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.text    = element_text(size = 7),
        plot.margin  = margin(4, 4, 4, 4)
      )
  })
  
  if (isTRUE(show_enrichment)) {
    enr <- meta_df$enrichment_score[match(motif_ids, meta_df$motif_id)]
    has_enr <- !all(is.na(enr))
    
    if (has_enr) {
      enr_df <- data.frame(
        motif_id = motif_ids,
        label    = vapply(motif_ids, function(m) meta_df$motif_name[meta_df$motif_id == m][1], character(1)),
        score    = enr,
        stringsAsFactors = FALSE
      )
      enr_df$label <- factor(enr_df$label, levels = rev(enr_df$label))
      
      bar_plot <- ggplot(enr_df, aes(x = score, y = label)) +
        geom_col(fill = "#378ADD", width = 0.6) +
        geom_text(aes(label = round(score, 2)), hjust = -0.1, size = 3) +
        xlab("Enrichment score") + ylab(NULL) +
        theme_minimal(base_size = 10) +
        theme(
          panel.grid.major.y = element_blank(),
          axis.text.y        = element_text(size = 9)
        )
      
      cowplot::plot_grid(
        cowplot::plot_grid(plotlist = logo_plots, ncol = ncol),
        bar_plot,
        ncol    = 2,
        rel_widths = c(3, 1)
      )
    } else {
      cowplot::plot_grid(plotlist = logo_plots, ncol = ncol)
    }
  } else {
    cowplot::plot_grid(plotlist = logo_plots, ncol = ncol)
  }
}

############################################### UI ###################################################

motif_plot_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
  tabPanel(
    HTML("Motif Plot"),
    h4("Transcription factor motif sequence logos"),
    "In this tab, users can browse and visualise sequence logos for transcription factor motifs ",
    "enriched in ATAC-seq peaks. Motifs can be filtered by name or enrichment score.",
    br(), br(),
    
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        
        h4("Filter motifs"),
        
        helper(
          textInput(
            ns("sc1mot_search"),
            "Search motif / TF name:",
            placeholder = "e.g. RUNX, MA0002"
          ),
          type    = "inline",
          content = "Type any part of the motif ID, motif name, or TF name. Case-insensitive."
        ),
        
        sliderInput(
          ns("sc1mot_enr_min"),
          "Minimum enrichment score:",
          min = 0, max = 10, value = 0, step = 0.1
        ),
        
        uiOutput(ns("sc1mot_enr_slider_note")),
        
        hr(),
        
        h4("Display"),
        
        uiOutput(ns("sc1mot_pick.ui")),
        
        sliderInput(
          ns("sc1mot_n"),
          "Number of motifs to show:",
          min = 1, max = 50, value = 10, step = 1
        ),
        
        sliderInput(
          ns("sc1mot_ncol"),
          "Columns:",
          min = 1, max = 8, value = 4, step = 1
        ),
        
        checkboxInput(ns("sc1mot_enr_bar"), "Show enrichment bar", value = TRUE),
        
        hr(),
        
        radioButtons(ns("sc1mot_psz"), "Plot size:",
                     choices = c("Small", "Medium", "Large"),
                     selected = "Medium", inline = TRUE)
      ),
      
      column(
        9,
        uiOutput(ns("sc1mot_oup.ui")),
        br(),
        downloadButton(ns("sc1mot_oup_pdf"), "Download PDF"),
        downloadButton(ns("sc1mot_oup_png"), "Download PNG"),
        br(),
        div(style = "display:inline-block",
            numericInput(ns("sc1mot_h"), "Height (in):", width = "120px", min = 2, max = 40, value = 10, step = 0.5)),
        div(style = "display:inline-block",
            numericInput(ns("sc1mot_w"), "Width (in):",  width = "120px", min = 2, max = 40, value = 14, step = 0.5)),
        hr(),
        h4("Motif table"),
        DTOutput(ns("sc1mot_table"))
      )
    )
  )
}

############################################### Server ###############################################

motif_plot_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    observe_helpers()
    
    motifs_rds  <- file.path(dir_inputs, "ATAC", "sc1motifs.rds")
    motifs_parq <- file.path(dir_inputs, "ATAC", "sc1motifs_meta.parquet")
    
    pwm_list <- if (file.exists(motifs_rds)) {
      readRDS(motifs_rds)
    } else {
      showNotification("sc1motifs.rds not found. Run prepShinyCellModular with ATAC assay.", type = "error")
      list()
    }
    
    meta_ds <- if (file.exists(motifs_parq)) arrow::open_dataset(motifs_parq) else NULL
    
    has_enrichment <- if (!is.null(meta_ds)) {
      enr <- meta_ds |> dplyr::select(enrichment_score) |> dplyr::collect()
      !all(is.na(enr$enrichment_score))
    } else {
      FALSE
    }
    
    output$sc1mot_enr_slider_note <- renderUI({
      if (!has_enrichment)
        tags$small(style = "color: #888;", "No enrichment scores available. Run FindMotifs and pass results to prepShinyCellModular.")
      else
        NULL
    })
    observe({
      if (has_enrichment) {
        max_enr <- meta_ds |>
          dplyr::summarise(m = max(enrichment_score, na.rm = TRUE)) |>
          dplyr::collect() |>
          dplyr::pull(m)
        updateSliderInput(session, "sc1mot_enr_min", max = round(max_enr, 1))
      }
    })
    filtered_meta <- reactive({
      req(!is.null(meta_ds))
      ds <- meta_ds
      if (has_enrichment && !is.null(input$sc1mot_enr_min) && input$sc1mot_enr_min > 0) {
        enr_min <- input$sc1mot_enr_min
        ds <- ds |> dplyr::filter(enrichment_score >= enr_min)
      }
      df <- as.data.frame(ds |> dplyr::collect())
      search <- trimws(input$sc1mot_search)
      if (nzchar(search)) {
        hits <- grepl(search, df$motif_id,   ignore.case = TRUE) |
          grepl(search, df$motif_name,  ignore.case = TRUE) |
          grepl(search, df$tf_name,     ignore.case = TRUE)
        df <- df[hits, , drop = FALSE]
      }
      df
    })
    output$sc1mot_pick.ui <- renderUI({
      df <- filtered_meta()
      choices <- setNames(df$motif_id, paste0(df$motif_name, " (", df$motif_id, ")"))
      selectizeInput(
        ns("sc1mot_pick"),
        "Select specific motifs (optional):",
        choices  = choices,
        selected = NULL,
        multiple = TRUE,
        options  = list(placeholder = "Leave empty to use top N", maxOptions = 200)
      )
    })
    active_ids <- reactive({
      df  <- filtered_meta()
      ids <- df$motif_id
      picked <- input$sc1mot_pick
      if (!is.null(picked) && length(picked) > 0) {
        ids <- intersect(picked, ids)
      } else {
        n <- min(input$sc1mot_n, length(ids))
        if (has_enrichment && !all(is.na(df$enrichment_score))) {
          df_ord <- df[order(df$enrichment_score, decreasing = TRUE, na.last = TRUE), ]
          ids    <- df_ord$motif_id[seq_len(n)]
        } else {
          ids <- ids[seq_len(n)]
        }
      }
      intersect(ids, names(pwm_list))
    })
    if (!exists("pList", inherits = TRUE)) {
      pList <<- c(Small = "350px", Medium = "550px", Large = "750px")
    }
    
    make_plot <- reactive({
      req(length(active_ids()) > 0)
      sc_motif_plot(
        pwm_list        = pwm_list,
        meta_df         = filtered_meta(),
        motif_ids       = active_ids(),
        ncol            = input$sc1mot_ncol,
        show_enrichment = isTRUE(input$sc1mot_enr_bar)
      )
    })
    
    output$sc1mot_oup <- renderPlot({ make_plot() })
    
    output$sc1mot_oup.ui <- renderUI({
      req(input$sc1mot_psz)
      plotOutput(ns("sc1mot_oup"), height = pList[input$sc1mot_psz])
    })
    
    output$sc1mot_table <- renderDT({
      df <- filtered_meta()
      ids <- active_ids()
      df <- df[df$motif_id %in% ids, , drop = FALSE]
      cols <- intersect(c("motif_id", "motif_name", "tf_name",
                          "enrichment_score", "pval", "padj"), colnames(df))
      df <- df[, cols, drop = FALSE]
      if ("enrichment_score" %in% colnames(df)) df$enrichment_score <- round(df$enrichment_score, 3)
      if ("pval"             %in% colnames(df)) df$pval             <- signif(df$pval, 3)
      if ("padj"             %in% colnames(df)) df$padj             <- signif(df$padj, 3)
      datatable(
        df,
        rownames  = FALSE,
        selection = "none",
        options   = list(pageLength = 10, scrollX = TRUE)
      )
    })
    
    output$sc1mot_oup_pdf <- downloadHandler(
      filename = function() paste0("motif_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content  = function(file) {
        ggsave(file, plot = make_plot(), device = "pdf",
               height = input$sc1mot_h, width = input$sc1mot_w, useDingbats = FALSE)
      }
    )
    
    output$sc1mot_oup_png <- downloadHandler(
      filename = function() paste0("motif_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) {
        ggsave(file, plot = make_plot(), device = "png",
               height = input$sc1mot_h, width = input$sc1mot_w)
      }
    )
  })
}

############################################### Registration ##########################################

register_tab(
  id     = "motif_plot",
  title  = "Motif Plot",
  ui     = motif_plot_ui,
  server = motif_plot_server
)