# id     = "proportions",
# title  = "Cells Proportions",

############################################### Functions ############################################

# Plot proportion plot
scProp <- function(inpConf, inpMeta,
                   inp1, inp2,
                   inpsub1, inpsub2,
                   inptyp, inpflp, inpfsz,
                   y_min = NULL, y_max = NULL,
                   x_order = NULL) {
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inp1]$ID,
    inpConf[UI == inp2]$ID,
    inpConf[UI == inpsub1]$ID
  ), with = FALSE]
  
  colnames(ggData) <- c("X", "grp", "sub")
  
  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2]
  }
  
  ggData <- ggData[, .(nCells = .N), by = c("X", "grp")]
  ggData <- ggData[, {
    tot <- sum(nCells)
    .SD[, .(
      pctCells = 100 * sum(nCells) / tot,
      nCells = nCells
    ), by = "grp"]
  }, by = "X"]
  
  # Factor grp and colors
  ggCol <- strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]]
  names(ggCol) <- levels(ggData$grp)
  ggLvl <- levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)]
  ggData$grp <- factor(ggData$grp, levels = ggLvl)
  ggCol <- ggCol[ggLvl]
  
  # Factor X order (independent from filter)
  if (!is.null(x_order) && length(x_order) > 0) {
    x_present <- unique(as.character(ggData$X))
    x_levels <- x_order[x_order %in% x_present]
    if (length(x_levels) > 0) {
      ggData$X <- factor(as.character(ggData$X), levels = x_levels)
    } else {
      ggData$X <- factor(ggData$X)
    }
  } else {
    ggData$X <- factor(ggData$X)
  }
  
  # Actual ggplot
  if (inptyp == "Proportion") {
    ggOut <- ggplot(ggData, aes(X, pctCells, fill = grp)) +
      geom_col() +
      ylab("Cell Proportion (%)")
  } else {
    ggOut <- ggplot(ggData, aes(X, nCells, fill = grp)) +
      geom_col() +
      ylab("Number of Cells")
  }
  
  if (isTRUE(inpflp)) {
    ggOut <- ggOut + coord_flip()
  }
  
  ggOut <- ggOut +
    xlab(inp1) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "right")
  
  # Fix Y scale (works with or without coord_flip)
  if (!is.null(y_min) && !is.null(y_max) && !is.na(y_min) && !is.na(y_max)) {
    ggOut <- ggOut + scale_y_continuous(limits = c(y_min, y_max))
  }
  
  ggOut
}

############################################## UI ################################################

scProp_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
  tabPanel(
    HTML("Proportion plot"),
    h4("Proportion / cell numbers across different cell information"),
    "In this tab, users can visualise the composition of single cells based on one discrete ",
    "cell information across another discrete cell information. ",
    "Usage examples include the library or cellcycle composition across clusters.",
    br(), br(),
    
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        
        selectInput(
          ns("sc1c2inp1"), "Cell information to plot (X-axis):",
          choices = sc1conf[grp == TRUE]$UI,
          selected = sc1def$grp2
        ),
        
        selectInput(
          ns("sc1c2inp2"), "Cell information to group / colour by:",
          choices = sc1conf[grp == TRUE]$UI,
          selected = sc1def$grp1
        ),
        
        radioButtons(
          ns("sc1c2typ"), "Plot value:",
          choices = c("Proportion", "CellNumbers"),
          selected = "Proportion", inline = TRUE
        ),
        
        checkboxInput(ns("sc1c2flp"), "Flip X/Y", value = FALSE),
        
        actionButton(ns("sc1c2togOrderX"), "Order X axis"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1c2togOrderX")),
          h5("Drag to reorder X axis groups"),
          uiOutput(ns("sc1c2xorder.ui")),
          actionButton(ns("sc1c2xorder_reset"), "Reset to default", class = "btn btn-primary")
        ),
        
        actionButton(ns("sc1c2fixscale"), "Fix Y scale", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("sc1c2fixscale")),
          fluidRow(
            column(6, numericInput(ns("sc1c2ymin"), "Y min", value = NULL, step = 0.1)),
            column(6, numericInput(ns("sc1c2ymax"), "Y max", value = NULL, step = 0.1))
          )
        ),
        
        actionButton(ns("sc1c2togL"), "Filter Cells"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1c2togL")),
          
          selectInput(
            ns("sc1c2sub1"), "Cell information to subset:",
            choices = sc1conf[grp == TRUE]$UI,
            selected = sc1def$grp1
          ),
          uiOutput(ns("sc1c2sub1.ui")),
          actionButton(ns("sc1c2sub1all"), "Select all groups", class = "btn btn-primary"),
          actionButton(ns("sc1c2sub1non"), "Deselect all groups", class = "btn btn-primary")
        ),
        
        br(), br(),
        
        actionButton(ns("sc1c2tog"), "Customize Plot"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1c2tog")),
          radioButtons(
            ns("sc1c2psz"), "Plot size:",
            choices = c("Small", "Medium", "Large"),
            selected = "Medium", inline = TRUE
          ),
          radioButtons(
            ns("sc1c2fsz"), "Font size:",
            choices = c("Small", "Medium", "Large"),
            selected = "Medium", inline = TRUE
          )
        )
      ),
      
      column(
        9,
        uiOutput(ns("sc1c2oup.ui")),
        downloadButton(ns("sc1c2oup.pdf"), "Download PDF"),
        downloadButton(ns("sc1c2oup.png"), "Download PNG"),
        br(),
        div(style = "display:inline-block",
            numericInput(ns("sc1c2oup.h"), "PDF / PNG height:", width = "138px",
                         min = 4, max = 20, value = 8, step = 0.5)),
        div(style = "display:inline-block",
            numericInput(ns("sc1c2oup.w"), "PDF / PNG width:", width = "138px",
                         min = 4, max = 20, value = 10, step = 0.5))
      )
    )
  )
}

############################################## Server ################################################

scProp_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    observe_helpers()
    optCrt <- "{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
    
    updateSelectizeInput(session, "sc1a1inp2", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene1, options = list(
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1a3inp1", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene1, options = list(
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1a3inp2", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene2, options = list(
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1b2inp1", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene1, options = list(
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1b2inp2", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene2, options = list(
                           maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1c1inp2", server = TRUE,
                         choices = c(sc1conf[is.na(fID)]$UI, names(sc1gene)),
                         selected = sc1conf[is.na(fID)]$UI[1], options = list(
                           maxOptions = length(sc1conf[is.na(fID)]$UI) + 3,
                           create = TRUE, persist = TRUE, render = I(optCrt)))
    
    if (!exists("pList2", inherits = TRUE)) {
      pList2 <<- c(Small = "450px", Medium = "650px", Large = "850px")
    }
    
    # Filter UI
    output$sc1c2sub1.ui <- renderUI({
      req(input$sc1c2sub1)
      sub <- strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(ns("sc1c2sub2"), "Select which cells to show", inline = TRUE,
                         choices = sub, selected = sub)
    })
    
    observeEvent(input$sc1c2sub1non, {
      req(input$sc1c2sub1)
      sub <- strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1c2sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = NULL, inline = TRUE)
    })
    
    observeEvent(input$sc1c2sub1all, {
      req(input$sc1c2sub1)
      sub <- strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1c2sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = sub, inline = TRUE)
    })
    
    # X axis order (independent from filter)
    x_levels_default <- reactive({
      req(input$sc1c2inp1)
      x_id <- sc1conf[UI == input$sc1c2inp1]$ID
      req(length(x_id) == 1, !is.na(x_id), x_id != "")
      levs <- levels(sc1meta[[x_id]])
      if (is.null(levs)) levs <- sort(unique(as.character(sc1meta[[x_id]])))
      levs
    })
    
    x_order <- reactiveVal(NULL)
    
    observeEvent(x_levels_default(), {
      x_order(x_levels_default())
    }, ignoreInit = FALSE)
    
    output$sc1c2xorder.ui <- renderUI({
      req(x_levels_default())
      ord <- x_order()
      if (is.null(ord) || length(ord) == 0) ord <- x_levels_default()
      ord <- ord[ord %in% x_levels_default()]
      if (length(ord) == 0) ord <- x_levels_default()
      
      sortable::rank_list(
        text = "X axis group order",
        labels = ord,
        input_id = ns("sc1c2xorder_rank")
      )
    })
    
    observeEvent(input$sc1c2xorder_rank, {
      req(input$sc1c2xorder_rank)
      x_order(input$sc1c2xorder_rank)
    })
    
    observeEvent(input$sc1c2xorder_reset, {
      x_order(x_levels_default())
    })
    
    x_order_final <- reactive({
      levs <- x_levels_default()
      ord <- x_order()
      if (is.null(ord) || length(ord) == 0) return(levs)
      ord <- ord[ord %in% levs]
      if (length(ord) == 0) levs else ord
    })
    
    # Plot output
    output$sc1c2oup <- renderPlot({
      scProp(
        sc1conf, sc1meta,
        input$sc1c2inp1, input$sc1c2inp2,
        input$sc1c2sub1, input$sc1c2sub2,
        input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz,
        y_min = input$sc1c2ymin,
        y_max = input$sc1c2ymax,
        x_order = x_order_final()
      )
    })
    
    output$sc1c2oup.ui <- renderUI({
      plotOutput(ns("sc1c2oup"), height = pList2[input$sc1c2psz])
    })
    
    # Downloads
    output$sc1c2oup.pdf <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1c2typ, "_", input$sc1c2inp1, "_", input$sc1c2inp2, ".pdf")
      },
      content = function(file) {
        ggsave(
          file, device = "pdf",
          height = input$sc1c2oup.h, width = input$sc1c2oup.w,
          useDingbats = FALSE,
          plot = scProp(
            sc1conf, sc1meta,
            input$sc1c2inp1, input$sc1c2inp2,
            input$sc1c2sub1, input$sc1c2sub2,
            input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz,
            y_min = input$sc1c2ymin,
            y_max = input$sc1c2ymax,
            x_order = x_order_final()
          )
        )
      }
    )
    
    output$sc1c2oup.png <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1c2typ, "_", input$sc1c2inp1, "_", input$sc1c2inp2, ".png")
      },
      content = function(file) {
        ggsave(
          file, device = "png",
          height = input$sc1c2oup.h, width = input$sc1c2oup.w,
          plot = scProp(
            sc1conf, sc1meta,
            input$sc1c2inp1, input$sc1c2inp2,
            input$sc1c2sub1, input$sc1c2sub2,
            input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz,
            y_min = input$sc1c2ymin,
            y_max = input$sc1c2ymax,
            x_order = x_order_final()
          )
        )
      }
    )
    
  })
}

############################################### Registration #################################################

register_tab(
  id     = "proportions",
  title  = "Cell Proportions",
  ui     = scProp_ui,
  server = scProp_server
)