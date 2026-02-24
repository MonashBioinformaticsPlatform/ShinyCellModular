# CellInfo3D vs GeneExpr3D
# This tab is base in the orginal "cellinfo_geneexpr" in ShinyCell
# id     = "cellinfo3D_geneexpr3D",
# title  = "CellInfo3D vs GeneExpr3D",

########################################## Functions ########################################################

# Return a data.table of available 3D reductions (prefixes) and their UI/ID mappings
scList3DReductions <- function(inpConf, inpMeta) {
  
  conf <- data.table::as.data.table(inpConf)
  meta_cols <- names(inpMeta)
  
  if (!("dimred" %in% names(conf))) stop("inpConf must contain a 'dimred' column.")
  if (!all(c("UI","ID") %in% names(conf))) stop("inpConf must contain columns UI and ID.")
  
  conf <- conf[dimred == TRUE]
  conf <- conf[!is.na(UI) & !is.na(ID)]
  conf <- conf[UI != "" & ID != ""]
  conf <- conf[ID %in% meta_cols]
  
  if (nrow(conf) == 0) {
    return(data.table::data.table())
  }
  
  conf[, dim := data.table::fifelse(grepl("([._-])?1$", UI), 1L,
                                    data.table::fifelse(grepl("([._-])?2$", UI), 2L,
                                                        data.table::fifelse(grepl("([._-])?3$", UI), 3L, NA_integer_)))]
  
  conf <- conf[!is.na(dim)]
  if (nrow(conf) == 0) {
    return(data.table::data.table())
  }
  
  conf[, base := gsub("([._-])?[123]$", "", UI)]
  
  bases <- conf[, .(
    has1 = any(dim == 1L),
    has2 = any(dim == 2L),
    has3 = any(dim == 3L)
  ), by = base][has1 & has2 & has3, .(base)]
  
  if (nrow(bases) == 0) {
    return(data.table::data.table())
  }
  
  out <- conf[base %in% bases$base, .(
    UI  = UI[1],
    ID  = ID[1]
  ), by = .(base, dim)]
  
  out_ui <- data.table::dcast(out, base ~ dim, value.var = "UI")
  out_id <- data.table::dcast(out, base ~ dim, value.var = "ID")
  
  data.table::setnames(out_ui, c("base", "UI_1", "UI_2", "UI_3"))
  data.table::setnames(out_id, c("base", "ID_1", "ID_2", "ID_3"))
  
  out_ui[out_id, on = "base"][]
}

scGet3DReduction <- function(inpConf, inpMeta, base = NULL) {
  
  dt <- scList3DReductions(inpConf, inpMeta)
  
  if (nrow(dt) == 0) {
    stop("No 3D reduction found among inpConf[dimred==TRUE] entries present in inpMeta.")
  }
  
  base_sel <- base
  if (is.null(base_sel) || !base_sel %in% dt$base) base_sel <- dt$base[1]
  
  row <- dt[base == base_sel][1]
  
  list(
    base = row$base,
    UI = list(X = row$UI_1, Y = row$UI_2, Z = row$UI_3),
    ID = list(X = row$ID_1, Y = row$ID_2, Z = row$ID_3)
  )
}

# Plot cell information on dimred 3D
scDRcell3D <- function(inpConf, inpMeta,
                       dr3d_base,
                       inp1, inpsub1, inpsub2,
                       inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab) {
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  red3 <- scGet3DReduction(inpConf, inpMeta, base = dr3d_base)
  
  ggData <- inpMeta[, c(
    red3$ID$X, red3$ID$Y, red3$ID$Z,
    inpConf[UI == inp1]$ID,
    inpConf[UI == inpsub1]$ID
  ), with = FALSE]
  
  data.table::setnames(ggData, c("X", "Y", "Z", "val", "sub"))
  
  bgCells <- FALSE
  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData  <- ggData[sub %in% inpsub2]
  }
  
  if (inpord == "Max-1st") {
    ggData <- ggData[order(val)]
  } else if (inpord == "Min-1st") {
    ggData <- ggData[order(-val)]
  } else if (inpord == "Random") {
    ggData <- ggData[sample(nrow(ggData))]
  }
  
  isDiscrete <- !is.na(inpConf[UI == inp1]$fCL)
  
  if (isDiscrete) {
    ggCol <- strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
    names(ggCol) <- levels(ggData$val)
    ggLvl <- levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)]
    ggData$val <- factor(ggData$val, levels = ggLvl)
    ggCol <- ggCol[ggLvl]
  }
  
  ggData$text <- as.character(ggData$val)
  if (bgCells) ggData2$text <- as.character(ggData2$val)
  
  hovertemplate <- paste(
    "Value: %{text}<br>",
    paste0(red3$UI$X, ": %{x:.2f}<br>"),
    paste0(red3$UI$Y, ": %{y:.2f}<br>"),
    paste0(red3$UI$Z, ": %{z:.2f}"),
    "<extra></extra>"
  )
  
  p <- plotly::plot_ly()
  
  if (bgCells) {
    p <- p %>%
      plotly::add_trace(
        data = ggData2,
        x = ~X, y = ~Y, z = ~Z,
        type = "scatter3d", mode = "markers",
        marker = list(size = inpsiz, color = "snow2"),
        text = ~text,
        hovertemplate = hovertemplate,
        showlegend = FALSE
      )
  }
  
  if (isDiscrete) {
    p <- p %>%
      plotly::add_trace(
        data = ggData,
        x = ~X, y = ~Y, z = ~Z,
        type = "scatter3d", mode = "markers",
        color = ~val,
        colors = unname(ggCol),
        text = ~text,
        marker = list(size = inpsiz),
        hovertemplate = hovertemplate
      )
  } else {
    p <- p %>%
      plotly::add_trace(
        data = ggData,
        x = ~X, y = ~Y, z = ~Z,
        type = "scatter3d", mode = "markers",
        text = ~text,
        marker = list(
          size = inpsiz,
          color = ggData$val,
          colorscale = cList[[inpcol]],
          showscale = TRUE
        ),
        hovertemplate = hovertemplate,
        showlegend = FALSE
      )
  }
  
  p <- p %>%
    plotly::layout(
      legend = list(itemsizing = "constant", font = list(size = 12)),
      scene = list(
        xaxis = list(title = red3$UI$X),
        yaxis = list(title = red3$UI$Y),
        zaxis = list(title = red3$UI$Z)
      )
    ) %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "cellinfo3D",
        width = 800,
        height = 600
      )
    )
  
  list(plot = p, reduction = red3$base)
}

# Plot gene expression on dimred 3D
scDRgene3D <- function(inpConf, inpMeta,
                       dr3d_base,
                       inp1, inpsub1, inpsub2,
                       inpH5, inpGene,
                       inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt) {
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  red3 <- scGet3DReduction(inpConf, inpMeta, base = dr3d_base)
  
  ggData <- inpMeta[, c(
    red3$ID$X, red3$ID$Y, red3$ID$Z,
    inpConf[UI == inpsub1]$ID
  ), with = FALSE]
  
  data.table::setnames(ggData, c("X", "Y", "Z", "sub"))
  
  h5file <- H5File$new(inpH5, mode = "r")
  on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
  h5data <- h5file[["grp"]][["data"]]
  
  ggData$val <- h5data$read(args = list(inpGene[inp1], quote(expr=)))
  ggData[val < 0]$val <- 0
  
  bgCells <- FALSE
  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData  <- ggData[sub %in% inpsub2]
  }
  
  if (inpord == "Max-1st") {
    ggData <- ggData[order(val)]
  } else if (inpord == "Min-1st") {
    ggData <- ggData[order(-val)]
  } else if (inpord == "Random") {
    ggData <- ggData[sample(nrow(ggData))]
  }
  
  ggData$text <- sprintf("%.3f", ggData$val)
  if (bgCells) ggData2$text <- sprintf("%.3f", ggData2$val)
  
  hovertemplate <- paste(
    paste0(inp1, ": %{text}<br>"),
    paste0(red3$UI$X, ": %{x:.2f}<br>"),
    paste0(red3$UI$Y, ": %{y:.2f}<br>"),
    paste0(red3$UI$Z, ": %{z:.2f}"),
    "<extra></extra>"
  )
  
  p <- plotly::plot_ly()
  
  if (bgCells) {
    p <- p %>%
      plotly::add_trace(
        data = ggData2,
        x = ~X, y = ~Y, z = ~Z,
        type = "scatter3d", mode = "markers",
        marker = list(size = inpsiz, color = "snow2"),
        text = ~text,
        hovertemplate = hovertemplate,
        showlegend = FALSE
      )
  }
  
  p <- p %>%
    plotly::add_trace(
      data = ggData,
      x = ~X, y = ~Y, z = ~Z,
      type = "scatter3d", mode = "markers",
      text = ~text,
      marker = list(
        size = inpsiz,
        color = ggData$val,
        colorscale = cList[[inpcol]],
        showscale = TRUE,
        colorbar = list(title = inp1)
      ),
      hovertemplate = hovertemplate,
      showlegend = FALSE
    ) %>%
    plotly::layout(
      legend = list(itemsizing = "constant", font = list(size = 12)),
      scene = list(
        xaxis = list(title = red3$UI$X),
        yaxis = list(title = red3$UI$Y),
        zaxis = list(title = red3$UI$Z)
      )
    ) %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "geneexpr3D",
        width = 800,
        height = 600
      )
    )
  
  list(plot = p, reduction = red3$base)
}


########################################## UI ########################################################

scDRnum3D_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
  tabPanel(
    title = HTML("CellInfo3D vs GeneExpr3D"),
    
    h4("Cell information vs gene expression on reduced dimensions 3D"),
    "In this tab, users can visualise both cell information and gene ",
    "expression side-by-side on low-dimensional representions in 3D.",
    br(), br(),
    
    fluidRow(
      column(
        3, h4("Dimension Reduction 3D"),
        selectInput(ns("sc1a1dr3d_base"), "reduction with more than 2 dimensions:", choices = character(0))
      ),
      
      column(
        3,
        actionButton(ns("sc1a1togL"), "Filter Cells"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a1togL")),
          selectInput(
            ns("sc1a1sub1"),
            "Cell information to subset:",
            choices = sc1conf[grp == TRUE]$UI,
            selected = sc1def$grp1
          ),
          uiOutput(ns("sc1a1sub1_ui")),
          actionButton(ns("sc1a1sub1all"), "Select all groups", class = "btn btn-primary"),
          actionButton(ns("sc1a1sub1non"), "Deselect all groups", class = "btn btn-primary")
        )
      ),
      
      column(
        6,
        actionButton(ns("sc1a1tog0"), "Customize plot"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a1tog0")),
          fluidRow(
            column(
              6,
              sliderInput(ns("sc1a1siz"), "Point size:", min = 0, max = 4, value = 1.25, step = 0.25),
              radioButtons(ns("sc1a1psz"), "Plot size:", choices = c("Small", "Medium", "Large"), selected = "Medium", inline = TRUE),
              radioButtons(ns("sc1a1fsz"), "Font size:", choices = c("Small", "Medium", "Large"), selected = "Medium", inline = TRUE)
            ),
            column(
              6,
              radioButtons(ns("sc1a1asp"), "Aspect ratio:", choices = c("Square", "Fixed", "Free"), selected = "Square", inline = TRUE),
              checkboxInput(ns("sc1a1txt"), "Show axis text", value = FALSE)
            )
          )
        )
      )
    ),
    
    fluidRow(
      column(
        6, style = "border-right: 2px solid black", h4("Cell information 3D"),
        fluidRow(
          column(
            6,
            selectInput(
              ns("sc1a1inp1"),
              "Cell information:",
              choices = sc1conf$UI,
              selected = sc1def$meta1
            )
          ),
          column(
            6,
            actionButton(ns("sc1a1tog1"), "Customize plot"),
            conditionalPanel(
              condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a1tog1")),
              radioButtons(ns("sc1a1col1"), "Colour (Continuous data):",
                           choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                           selected = "Blue-Yellow-Red"),
              radioButtons(ns("sc1a1ord1"), "Plot order:",
                           choices = c("Max-1st", "Min-1st", "Original", "Random"),
                           selected = "Original", inline = TRUE),
              checkboxInput(ns("sc1a1lab1"), "Show cell info labels", value = TRUE)
            )
          )
        ),
        fluidRow(column(12, uiOutput(ns("sc1a1oup1_ui"))))
      ),
      
      column(
        6, h4("Gene expression 3D"),
        fluidRow(
          column(
            6,
            selectInput(ns("sc1a1inp2"), "Gene name:", choices = NULL)
          ),
          column(
            6,
            actionButton(ns("sc1a1tog2"), "Customize plot"),
            conditionalPanel(
              condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a1tog2")),
              radioButtons(ns("sc1a1col2"), "Colour:",
                           choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                           selected = "White-Red"),
              radioButtons(ns("sc1a1ord2"), "Plot order:",
                           choices = c("Max-1st", "Min-1st", "Original", "Random"),
                           selected = "Max-1st", inline = TRUE)
            )
          )
        ),
        fluidRow(column(12, uiOutput(ns("sc1a1oup2_ui"))))
      )
    )
  )
}

########################################## Server ########################################################

scDRnum3D_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    observe_helpers()
    
    optCrt <- "{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
    updateSelectizeInput(
      session, "sc1a1inp2",
      choices = names(sc1gene),
      server = TRUE,
      selected = sc1def$gene1,
      options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))
    )
    
    if (!exists("pList", inherits = TRUE)) {
      pList <<- c(Small = "350px", Medium = "550px", Large = "750px")
    }
    
    dr3d_tbl <- reactive({
      scList3DReductions(sc1conf, sc1meta)
    })
    
    observeEvent(dr3d_tbl(), {
      tbl <- dr3d_tbl()
      choices <- tbl$base
      if (length(choices) == 0) choices <- character(0)
      
      current <- isolate(input$sc1a1dr3d_base)
      selected <- if (!is.null(current) && current %in% choices) current else if (length(choices) > 0) choices[1] else character(0)
      
      updateSelectInput(session, "sc1a1dr3d_base", choices = choices, selected = selected)
    }, ignoreInit = FALSE)
    
    output$sc1a1sub1_ui <- renderUI({
      req(input$sc1a1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(ns("sc1a1sub2"), "Select which cells to show", inline = TRUE,
                         choices = sub, selected = sub)
    })
    
    observeEvent(input$sc1a1sub1non, {
      req(input$sc1a1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1a1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = NULL, inline = TRUE)
    })
    
    observeEvent(input$sc1a1sub1all, {
      req(input$sc1a1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1a1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = sub, inline = TRUE)
    })
    
    plot1 <- reactive({
      req(input$sc1a1dr3d_base, input$sc1a1inp1)
      scDRcell3D(
        sc1conf, sc1meta,
        dr3d_base = input$sc1a1dr3d_base,
        inp1    = input$sc1a1inp1,
        inpsub1 = input$sc1a1sub1,
        inpsub2 = input$sc1a1sub2,
        inpsiz  = input$sc1a1siz,
        inpcol  = input$sc1a1col1,
        inpord  = input$sc1a1ord1,
        inpfsz  = input$sc1a1fsz,
        inpasp  = input$sc1a1asp,
        inptxt  = input$sc1a1txt,
        inplab  = input$sc1a1lab1
      )
    })
    
    plot2 <- reactive({
      req(input$sc1a1dr3d_base, input$sc1a1inp2)
      scDRgene3D(
        sc1conf, sc1meta,
        dr3d_base = input$sc1a1dr3d_base,
        inp1    = input$sc1a1inp2,
        inpsub1 = input$sc1a1sub1,
        inpsub2 = input$sc1a1sub2,
        inpH5   = file.path(dir_inputs, "sc1gexpr.h5"),
        inpGene = sc1gene,
        inpsiz  = input$sc1a1siz,
        inpcol  = input$sc1a1col2,
        inpord  = input$sc1a1ord2,
        inpfsz  = input$sc1a1fsz,
        inpasp  = input$sc1a1asp,
        inptxt  = input$sc1a1txt
      )
    })
    
    output$sc1a1oup1 <- plotly::renderPlotly({
      plot1()$plot
    })
    
    output$sc1a1oup1_ui <- renderUI({
      req(input$sc1a1psz)
      plotly::plotlyOutput(ns("sc1a1oup1"), height = pList[input$sc1a1psz])
    })
    
    output$sc1a1oup2 <- plotly::renderPlotly({
      plot2()$plot
    })
    
    output$sc1a1oup2_ui <- renderUI({
      req(input$sc1a1psz)
      plotly::plotlyOutput(ns("sc1a1oup2"), height = pList[input$sc1a1psz])
    })
    
  })
}

############################################### Registration #################################################

register_tab(
  id     = "cellinfo3D_geneexpr3D",
  title  = "CellInfo3D vs GeneExpr3D",
  ui     = scDRnum3D_ui,
  server = scDRnum3D_server
)