# id     = "violinbox",
# title  = "Violin / BoxPlot",

######################################################## Functions ######################################################

scVioBox <- function(inpConf, inpMeta, inp1, inp2,
                     inpsub1, inpsub2, inpH5, inpGene,
                     inptyp, inppts, inpsiz, inpfsz, inscale_min, inscale_max,
                     x_order = NULL) {
  
  expr_min <- inscale_min
  expr_max <- inscale_max
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  ggData <- inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), with = FALSE]
  colnames(ggData) <- c("X", "sub")
  
  if (inp2 %in% inpConf$UI) {
    ggData$val <- inpMeta[[inpConf[UI == inp2]$ID]]
  } else {
    h5file <- H5File$new(inpH5, mode = "r")
    on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
    h5data <- h5file[["grp"]][["data"]]
    ggData$val <- h5data$read(args = list(inpGene[inp2], quote(expr=)))
    ggData[val < 0]$val <- 0
    set.seed(42)
    tmpNoise <- rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000
    ggData$val <- ggData$val + tmpNoise
  }
  
  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2]
  }
  
  ggCol <- strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
  names(ggCol) <- levels(ggData$X)
  
  ggLvl <- levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)]
  
  if (!is.null(x_order) && length(x_order) > 0) {
    ggLvl <- x_order[x_order %in% ggLvl]
  }
  
  ggData$X <- factor(ggData$X, levels = ggLvl)
  ggCol <- ggCol[ggLvl]
  
  if (inptyp == "violin") {
    ggOut <- ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width")
  } else {
    ggOut <- ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot()
  }
  
  if (isTRUE(inppts)) {
    ggOut <- ggOut + geom_jitter(size = inpsiz, shape = 16)
  }
  
  ggOut <- ggOut +
    xlab(inp1) + ylab(inp2) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
  
  if (!is.null(expr_min) && !is.null(expr_max) &&
      !is.na(expr_min) && !is.na(expr_max)) {
    ggOut <- ggOut + scale_y_continuous(limits = c(expr_min, expr_max))
  }
  
  ggOut
}

scVioSummary <- function(inpConf, inpMeta, inp1, inp2,
                         inpsub1, inpsub2, inpH5, inpGene, inpsplt,
                         x_order = NULL) {
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  ggData <- inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), with = FALSE]
  colnames(ggData) <- c("group", "sub")
  
  is_meta_y <- inp2 %in% inpConf$UI
  
  if (is_meta_y) {
    ggData$val <- inpMeta[[ inpConf[UI == inp2]$ID ]]
  } else {
    if (is.null(inpGene) || !(inp2 %in% names(inpGene))) {
      stop(paste0("Selected Y (", inp2, ") is not a meta variable and is not in sc1gene."))
    }
    h5file <- H5File$new(inpH5, mode = "r")
    on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
    h5data <- h5file[["grp"]][["data"]]
    ggData$val <- h5data$read(args = list(inpGene[inp2], quote(expr=)))
    ggData[val < 0]$val <- 0
  }
  
  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2]
  }
  
  if (is.na(inpConf[UI == inp1]$fCL)) {
    if (inpsplt == "Quartile") nBk <- 4
    if (inpsplt == "Decile")   nBk <- 10
    ggData$group <- cut(ggData$group, breaks = nBk)
  }
  
  ggData$express <- ggData$val > 0
  ggData$express[is.na(ggData$express)] <- FALSE
  
  summaryTable <- ggData[, .(
    nCells     = .N,
    nExpress   = sum(express, na.rm = TRUE),
    pctExpress = 100 * sum(express, na.rm = TRUE) / .N,
    minVal     = suppressWarnings(min(val, na.rm = TRUE)),
    maxVal     = suppressWarnings(max(val, na.rm = TRUE)),
    meanVal    = suppressWarnings(mean(val, na.rm = TRUE))
  ), by = "group"]
  
  if (!is.null(x_order) && length(x_order) > 0) {
    summaryTable$group <- factor(summaryTable$group, levels = x_order)
    summaryTable <- summaryTable[order(summaryTable$group)]
  } else {
    summaryTable <- summaryTable[order(group)]
  }
  
  colnames(summaryTable)[which(names(summaryTable) == "minVal")]  <- paste0("min_", inp2)
  colnames(summaryTable)[which(names(summaryTable) == "maxVal")]  <- paste0("max_", inp2)
  colnames(summaryTable)[which(names(summaryTable) == "meanVal")] <- paste0("mean_", inp2)
  
  summaryTable[]
}

######################################################## UI ######################################################

scVioBox_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
  tabPanel(
    HTML("Violinplot / Boxplot"),
    h4("Cell information / gene expression violin plot / box plot"),
    "In this tab, users can visualise the gene expression or continuous cell information ",
    "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).",
    br(), br(),
    
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        
        selectInput(
          ns("sc1c1inp1"),
          "Cell information (X-axis):",
          choices = sc1conf[grp == TRUE]$UI,
          selected = sc1def$grp1
        ),
        
        selectInput(ns("sc1c1inp2"), "Cell Info / Gene name (Y-axis):", choices = NULL),
        
        radioButtons(ns("sc1c1typ"), "Plot type:",
                     choices = c("violin", "boxplot"),
                     selected = "violin", inline = TRUE),
        
        checkboxInput(ns("sc1c1pts"), "Show data points", value = FALSE),
        
        actionButton(ns("sc1c1togOrderX"), "Order X axis"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1c1togOrderX")),
          h5("Drag to reorder X axis groups"),
          uiOutput(ns("sc1c1xorder.ui")),
          actionButton(ns("sc1c1xorder_reset"), "Reset to default", class = "btn btn-primary")
        ),
        
        br(), br(),
        
        actionButton(ns("sc1c1togL"), "Filter Cells"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1c1togL")),
          selectInput(ns("sc1c1sub1"), "Cell information to subset:",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp1),
          uiOutput(ns("sc1c1sub1.ui")),
          actionButton(ns("sc1c1sub1all"), "Select all groups", class = "btn btn-primary"),
          actionButton(ns("sc1c1sub1non"), "Deselect all groups", class = "btn btn-primary")
        ),
        
        br(), br(),
        
        actionButton(ns("sc1c1tog"), "Customize Plot"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1c1tog")),
          sliderInput(ns("sc1c1siz"), "Data point size:",
                      min = 0, max = 4, value = 1.25, step = 0.25),
          radioButtons(ns("sc1c1psz"), "Plot size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Medium", inline = TRUE),
          radioButtons(ns("sc1c1fsz"), "Font size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Medium", inline = TRUE)
        ),
        
        actionButton(ns("sc1c1fixscale"), "Fix Y scale", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("sc1c1fixscale")),
          fluidRow(
            column(6, numericInput(ns("sc1c1ymin"), "Y min", value = NULL, step = 0.1)),
            column(6, numericInput(ns("sc1c1ymax"), "Y max", value = NULL, step = 0.1))
          )
        )
      ),
      
      column(
        9,
        uiOutput(ns("sc1c1oup.ui")),
        downloadButton(ns("sc1c1oup.pdf"), "Download PDF"),
        downloadButton(ns("sc1c1oup.png"), "Download PNG"),
        br(),
        
        div(style = "display:inline-block",
            numericInput(ns("sc1c1oup.h"), "PDF / PNG height:", width = "138px",
                         min = 4, max = 20, value = 8, step = 0.5)),
        
        div(style = "display:inline-block",
            numericInput(ns("sc1c1oup.w"), "PDF / PNG width:", width = "138px",
                         min = 4, max = 20, value = 10, step = 0.5)),
        
        br(),
        
        actionButton(ns("sc1c1tog1"), "Show cell numbers / statistics"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1c1tog1")),
          h4("Cell numbers / statistics"),
          radioButtons(ns("sc1c1splt"), "Split continuous cell info into:",
                       choices = c("Quartile", "Decile"),
                       selected = "Decile", inline = TRUE),
          DT::DTOutput(ns("sc1c1.dt"))
        )
      )
    )
  )
}

######################################################## Server ######################################################

scVioBox_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    observe_helpers()
    
    optCrt <- "{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
    
    updateSelectizeInput(
      session, "sc1c1inp2", server = TRUE,
      choices = c(sc1conf[is.na(fID)]$UI, names(sc1gene)),
      selected = sc1conf[is.na(fID)]$UI[1],
      options = list(
        maxOptions = length(sc1conf[is.na(fID)]$UI) + 3,
        create = TRUE, persist = TRUE, render = I(optCrt)
      )
    )
    
    if (!exists("pList2", inherits = TRUE)) {
      pList2 <<- c(Small = "450px", Medium = "650px", Large = "850px")
    }
    
    # Filter UI (independent)
    output$sc1c1sub1.ui <- renderUI({
      req(input$sc1c1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(
        ns("sc1c1sub2"), "Select which cells to show",
        inline = TRUE,
        choices = sub,
        selected = sub
      )
    })
    
    observeEvent(input$sc1c1sub1non, {
      req(input$sc1c1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(
        session,
        inputId = "sc1c1sub2",
        label = "Select which cells to show",
        choices = sub,
        selected = NULL,
        inline = TRUE
      )
    })
    
    observeEvent(input$sc1c1sub1all, {
      req(input$sc1c1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(
        session,
        inputId = "sc1c1sub2",
        label = "Select which cells to show",
        choices = sub,
        selected = sub,
        inline = TRUE
      )
    })
    
    # X axis ordering levels come from the X axis variable itself (independent from filter)
    x_levels_default <- reactive({
      req(input$sc1c1inp1)
      x_id <- sc1conf[UI == input$sc1c1inp1]$ID
      req(length(x_id) == 1, !is.na(x_id), x_id != "")
      levs <- levels(sc1meta[[x_id]])
      if (is.null(levs)) levs <- sort(unique(as.character(sc1meta[[x_id]])))
      levs
    })
    
    x_order <- reactiveVal(NULL)
    
    observeEvent(x_levels_default(), {
      x_order(x_levels_default())
    }, ignoreInit = FALSE)
    
    output$sc1c1xorder.ui <- renderUI({
      req(x_levels_default())
      ord <- x_order()
      if (is.null(ord) || length(ord) == 0) ord <- x_levels_default()
      ord <- ord[ord %in% x_levels_default()]
      if (length(ord) == 0) ord <- x_levels_default()
      
      sortable::rank_list(
        text = "X axis group order",
        labels = ord,
        input_id = ns("sc1c1xorder_rank")
      )
    })
    
    observeEvent(input$sc1c1xorder_rank, {
      req(input$sc1c1xorder_rank)
      x_order(input$sc1c1xorder_rank)
    })
    
    observeEvent(input$sc1c1xorder_reset, {
      x_order(x_levels_default())
    })
    
    x_order_final <- reactive({
      levs <- x_levels_default()
      ord <- x_order()
      if (is.null(ord) || length(ord) == 0) return(levs)
      ord <- ord[ord %in% levs]
      if (length(ord) == 0) levs else ord
    })
    
    output$sc1c1oup <- renderPlot({
      req(input$sc1c1inp1, input$sc1c1inp2)
      scVioBox(
        sc1conf, sc1meta,
        input$sc1c1inp1, input$sc1c1inp2,
        input$sc1c1sub1, input$sc1c1sub2,
        file.path(dir_inputs, "sc1gexpr.h5"),
        sc1gene,
        input$sc1c1typ, input$sc1c1pts,
        input$sc1c1siz, input$sc1c1fsz,
        input$sc1c1ymin, input$sc1c1ymax,
        x_order = x_order_final()
      )
    })
    
    output$sc1c1oup.ui <- renderUI({
      req(input$sc1c1psz)
      plotOutput(ns("sc1c1oup"), height = pList2[input$sc1c1psz])
    })
    
    output$sc1c1.dt <- DT::renderDT({
      req(input$sc1c1inp1, input$sc1c1inp2)
      
      ggData <- scVioSummary(
        sc1conf, sc1meta,
        input$sc1c1inp1, input$sc1c1inp2,
        input$sc1c1sub1, input$sc1c1sub2,
        file.path(dir_inputs, "sc1gexpr.h5"),
        sc1gene,
        input$sc1c1splt,
        x_order = x_order_final()
      )
      
      DT::datatable(
        ggData,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))
      ) %>%
        DT::formatRound(columns = "pctExpress", digits = 2)
    }, server = FALSE)
    
    output$sc1c1oup.pdf <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1c1typ, "_", input$sc1c1inp1, "_", input$sc1c1inp2, ".pdf")
      },
      content = function(file) {
        ggsave(
          file,
          device = "pdf",
          height = input$sc1c1oup.h,
          width  = input$sc1c1oup.w,
          useDingbats = FALSE,
          plot = scVioBox(
            sc1conf, sc1meta,
            input$sc1c1inp1, input$sc1c1inp2,
            input$sc1c1sub1, input$sc1c1sub2,
            file.path(dir_inputs, "sc1gexpr.h5"),
            sc1gene,
            input$sc1c1typ, input$sc1c1pts,
            input$sc1c1siz, input$sc1c1fsz,
            input$sc1c1ymin, input$sc1c1ymax,
            x_order = x_order_final()
          )
        )
      }
    )
    
    output$sc1c1oup.png <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1c1typ, "_", input$sc1c1inp1, "_", input$sc1c1inp2, ".png")
      },
      content = function(file) {
        ggsave(
          file,
          device = "png",
          height = input$sc1c1oup.h,
          width  = input$sc1c1oup.w,
          plot = scVioBox(
            sc1conf, sc1meta,
            input$sc1c1inp1, input$sc1c1inp2,
            input$sc1c1sub1, input$sc1c1sub2,
            file.path(dir_inputs, "sc1gexpr.h5"),
            sc1gene,
            input$sc1c1typ, input$sc1c1pts,
            input$sc1c1siz, input$sc1c1fsz,
            input$sc1c1ymin, input$sc1c1ymax,
            x_order = x_order_final()
          )
        )
      }
    )
  })
}

############################################### Registration #################################################

register_tab(
  id     = "violin_boxplot",
  title  = "Violin / BoxPlot",
  ui     = scVioBox_ui,
  server = scVioBox_server
)