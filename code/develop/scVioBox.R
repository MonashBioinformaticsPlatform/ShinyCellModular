######################################################## Function ######################################################

# Plot violin / boxplot 
scVioBox <- function(inpConf, inpMeta, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inptyp, inppts, inpsiz, inpfsz,inscale_min,inscale_max){ 
  expr_min <- inscale_min               # or a value you choose
  expr_max <- inscale_max
  
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("X", "sub") 
  
  # Load in either cell meta or gene expr
  if(inp2 %in% inpConf$UI){ 
    ggData$val = inpMeta[[inpConf[UI == inp2]$ID]] 
  } else { 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
    ggData[val < 0]$val = 0 
    set.seed(42) 
    tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000 
    ggData$val = ggData$val + tmpNoise 
    h5file$close_all() 
  } 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$X) 
  ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)] 
  ggData$X = factor(ggData$X, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "violin"){ 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width") 
  } else { 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot() 
  } 
  if(inppts){ 
    ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16) 
  } 
  
  ggOut = ggOut + xlab(inp1) + ylab(inp2) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
  
  if (!is.null(expr_min) && !is.null(expr_max) &&
      !is.na(expr_min) && !is.na(expr_max)) {
    ggOut <- ggOut + scale_y_continuous(limits = c(expr_min, expr_max))
  }
  
  return(ggOut) 
} 
 
scVioSummary <- function(inpConf, inpMeta, inp1, inp2, 
                         inpsub1, inpsub2, inpH5, inpGene, inpsplt){ 
  if(is.null(inpsub1)){ inpsub1 = inpConf$UI[1] } 
  
  # Prepare ggData
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  
  colnames(ggData) = c("group", "sub") 
  
  # Load expression values
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  h5file$close_all() 
  
  # Set negative values to zero
  ggData[val < 0]$val = 0 
  
  # Filter subset if necessary
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Split group if necessary
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    if(inpsplt == "Quartile"){ nBk = 4 } 
    if(inpsplt == "Decile"){ nBk = 10 } 
    ggData$group = cut(ggData$group, breaks = nBk) 
  } 
  
  # Compute statistics per group
  ggData$express = ggData$val > 0 
  
  summaryTable = ggData[, .(
    nCells = .N,
    nExpress = sum(express),
    minExpr = ifelse(sum(express) > 0, min(val[express == TRUE]), 0),
    maxExpr = ifelse(sum(express) > 0, max(val[express == TRUE]), 0),
    meanExpr = ifelse(sum(express) > 0, mean(val[express == TRUE]), 0)
  ), by = "group"]
  
  # Compute % expressing cells
  summaryTable[, pctExpress := 100 * nExpress / nCells]
  
  # Order by group
  summaryTable = summaryTable[order(group)]
  
  # Rename mean expression column for clarity
  colnames(summaryTable)[5] = paste0("meanExpr_", inp2) 
  
  return(summaryTable)
  #return( inp1)
}

######################################################## UI ######################################################

scVioBox_ui <- function(id, sc1conf, sc1def) {
  tabPanel( 
      HTML("Violinplot / Boxplot"),  
    h4("Cell information / gene expression violin plot / box plot"), 
    "In this tab, users can visualise the gene expression or continuous cell information ",  
    "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        selectInput("sc1c1inp1", "Cell information (X-axis):", 
                    choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1def$grp1) %>%  
          helper(type = "inline", size = "m", fade = TRUE, 
                  title = "Cell information to group cells by",  
                  content = c("Select categorical cell information to group cells by",  
                              "- Single cells are grouped by this categorical covariate",  
                              "- Plotted as the X-axis of the violin plot / box plot")),  
        selectInput("sc1c1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
          helper(type = "inline", size = "m", fade = TRUE, 
                  title = "Cell Info / Gene to plot", 
                  content = c("Select cell info / gene to plot on Y-axis", 
                              "- Can be continuous cell information (e.g. nUMIs / scores)", 
                              "- Can also be gene expression")), 
        radioButtons("sc1c1typ", "Plot type:", 
                      choices = c("violin", "boxplot"), 
                      selected = "violin", inline = TRUE), 
        checkboxInput("sc1c1pts", "Show data points", value = FALSE), 
        actionButton("sc1c1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1c1togL % 2 == 1", 
          selectInput("sc1c1sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1c1sub1.ui"), 
          actionButton("sc1c1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1c1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("sc1c1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1c1tog % 2 == 1", 
          sliderInput("sc1c1siz", "Data point size:",  
                      min = 0, max = 4, value = 1.25, step = 0.25),  
          radioButtons("sc1c1psz", "Plot size:", 
                        choices = c("Small", "Medium", "Large"), 
                        selected = "Medium", inline = TRUE), 
          radioButtons("sc1c1fsz", "Font size:", 
                        choices = c("Small", "Medium", "Large"), 
                        selected = "Medium", inline = TRUE)) 
          checkboxInput("sc1c1fixscale", "Fix Y scale", value = FALSE),
          conditionalPanel(
            condition = "input.sc1c1fixscale == true",
            fluidRow(
              column(
                6,
                numericInput("sc1c1ymin", "Y min", value = NULL, step = 0.1)
              ),
              column(
                6,
                numericInput("sc1c1ymax", "Y max", value = NULL, step = 0.1)
              )
            )
          )
      ), # End of column (6 space) 
      column(9, uiOutput("sc1c1oup.ui"),  
              downloadButton("sc1c1oup.pdf", "Download PDF"),  
              downloadButton("sc1c1oup.png", "Download PNG"), br(), 
              div(style="display:inline-block", 
                  numericInput("sc1c1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 8, step = 0.5)), 
              div(style="display:inline-block", 
                  numericInput("sc1c1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)),
              ####################################################################
              br(), 
              actionButton("sc1c1tog1", "Toggle to show cell numbers / statistics"), 
            conditionalPanel( 
                condition = "input.sc1c1tog1 % 2 == 1", 
                h4("Cell numbers / statistics"), 
                radioButtons("sc1c1splt", "Split continuous cell info into:", 
                            choices = c("Quartile", "Decile"), 
                            selected = "Decile", inline = TRUE), 
                dataTableOutput("sc1c1.dt")
              ####################################################################
        )  # End of column (6 space) 
      ) 
    )    # End of fluidRow (4 space) 
  )     # End of tab (2 space) 
 }

######################################################## Server ######################################################
scVioBox_server <- function(id, sc1conf, sc1def, sc1meta, sc1gene, dir_inputs) {

  ## UI for subset selection
  output$sc1c1sub1.ui <- renderUI({
    sub <- strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]]
    checkboxGroupInput(
      "sc1c1sub2", "Select which cells to show",
      inline  = TRUE,
      choices = sub,
      selected = sub
    )
  })

  observeEvent(input$sc1c1sub1non, {
    sub <- strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1c1sub2",
      label   = "Select which cells to show",
      choices = sub,
      selected = NULL,
      inline  = TRUE
    )
  })

  observeEvent(input$sc1c1sub1all, {
    sub <- strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1c1sub2",
      label   = "Select which cells to show",
      choices = sub,
      selected = sub,
      inline  = TRUE
    )
  })

  ## Main plot
  output$sc1c1oup <- renderPlot({
    scVioBox(
      sc1conf, sc1meta,
      input$sc1c1inp1, input$sc1c1inp2,
      input$sc1c1sub1, input$sc1c1sub2,
      paste0(dir_inputs, "sc1gexpr.h5"),
      sc1gene,
      input$sc1c1typ, input$sc1c1pts,
      input$sc1c1siz, input$sc1c1fsz,
      input$sc1c1ymin, input$sc1c1ymax
    )
  })

  output$sc1c1oup.ui <- renderUI({
    plotOutput("sc1c1oup", height = pList2[input$sc1c1psz])
  })

  ## Summary table
  output$sc1c1.dt <- renderDataTable({
    # optional debug
    # print(input$sc1c1inp1)
    # print(input$sc1c1inp2)

    ggData <- scVioSummary(
      sc1conf, sc1meta,
      input$sc1c1inp1, input$sc1c1inp2,
      input$sc1c1sub1, input$sc1c1sub2,
      paste0(dir_inputs, "sc1gexpr.h5"),
      sc1gene,
      input$sc1c1splt
    )

    datatable(
      ggData,
      rownames   = FALSE,
      extensions = "Buttons",
      options    = list(
        pageLength = -1,
        dom        = "tB",
        buttons    = c("copy", "csv", "excel")
      )
    ) %>%
      formatRound(columns = c("pctExpress"), digits = 2)
  })

  ## PDF download
  output$sc1c1oup.pdf <- downloadHandler(
    filename = function() {
      paste0("sc1", input$sc1c1typ, "_", input$sc1c1inp1, "_", input$sc1c1inp2, ".pdf")
    },
    content = function(file) {
      ggsave(
        file,
        device     = "pdf",
        height     = input$sc1c1oup.h,
        width      = input$sc1c1oup.w,
        useDingbats = FALSE,
        plot = scVioBox(
          sc1conf, sc1meta,
          input$sc1c1inp1, input$sc1c1inp2,
          input$sc1c1sub1, input$sc1c1sub2,
          paste0(dir_inputs, "sc1gexpr.h5"),
          sc1gene,
          input$sc1c1typ, input$sc1c1pts,
          input$sc1c1siz, input$sc1c1fsz,
          input$sc1c1ymin, input$sc1c1ymax
        )
      )
    }
  )

  ## PNG download
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
          paste0(dir_inputs, "sc1gexpr.h5"),
          sc1gene,
          input$sc1c1typ, input$sc1c1pts,
          input$sc1c1siz, input$sc1c1fsz,
          input$sc1c1ymin, input$sc1c1ymax
        )
      )
    }
  )
}
