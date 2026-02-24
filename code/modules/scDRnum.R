# This tab is base in the orginal "cellinfo_geneexpr" in ShinyCell
# id     = "cellinfo_geneexpr",
# title  = "CellInfo vs GeneExpr",

############################################### Functions ############################################


scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2,
                    inpH5, inpGene, inpsplt){
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  ggData <- inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), with = FALSE]
  colnames(ggData) <- c("group", "sub")
  
  h5file <- H5File$new(inpH5, mode = "r")
  on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
  h5data <- h5file[["grp"]][["data"]]
  
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr=)))
  ggData[val2 < 0]$val2 <- 0
  
  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2]
  }
  
  if (is.na(inpConf[UI == inp1]$fCL)) {
    if (inpsplt == "Quartile") nBk <- 4
    if (inpsplt == "Decile")   nBk <- 10
    ggData$group <- cut(ggData$group, breaks = nBk)
  }
  
  ggData$express <- FALSE
  ggData[val2 > 0]$express <- TRUE
  ggData1 <- ggData[express == TRUE, .(nExpress = .N), by = "group"]
  ggData  <- ggData[, .(nCells = .N), by = "group"]
  ggData  <- ggData1[ggData, on = "group"]
  ggData  <- ggData[, c("group", "nCells", "nExpress"), with = FALSE]
  ggData[is.na(nExpress)]$nExpress <- 0
  ggData$pctExpress <- 100 * ggData$nExpress / ggData$nCells
  ggData <- ggData[order(group)]
  colnames(ggData)[3] <- paste0(colnames(ggData)[3], "_", inp2)
  ggData
}



# Plot gene expression on dimred 
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val < 0]$val = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +  
    scale_color_gradientn(inp1, colours = cList[[inpcol]]) + 
    guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 

# Plot cell information on dimred 
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "val", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Do factoring if required 
  if(!is.na(inpConf[UI == inp1]$fCL)){ 
    ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
    names(ggCol) = levels(ggData$val) 
    ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
    ggData$val = factor(ggData$val, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    ggOut = ggOut + scale_color_gradientn("", colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  } else { 
    sListX = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
    sListX = 0.75 * (sList - (1.5 * floor(sListX/50))) 
    ggOut = ggOut + scale_color_manual("", values = ggCol) + 
      guides(color = guide_legend(override.aes = list(size = 5),  
                                  nrow = inpConf[UI == inp1]$fRow)) + 
      theme(legend.text = element_text(size = sListX[inpfsz])) 
    if(inplab){ 
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
      lListX = min(nchar(paste0(ggData3$val, collapse = "")), 200) 
      lListX = lList - (0.25 * floor(lListX/50)) 
      ggOut = ggOut + 
        geom_text_repel(data = ggData3, aes(X, Y, label = val), 
                        color = "grey10", bg.color = "grey95", bg.r = 0.15, 
                        size = lListX[inpfsz], seed = 42) 
    } 
  } 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
}

############################################### UI #####################################################

scDRnum_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
  tabPanel(
    title = HTML("CellInfo vs GeneExpr"),
    
    h4("Cell information vs gene expression on reduced dimensions"),
    "In this tab, users can visualise both cell information and gene ",
    "expression side-by-side on low-dimensional representions.",
    br(), br(),
    
    fluidRow(
      column(
        3, h4("Dimension Reduction"),
        fluidRow(
          column(
            12,
            selectInput(ns("sc1a1drX"), "X-axis:", choices = sc1conf[dimred == TRUE]$UI, selected = sc1def$dimred[1]),
            selectInput(ns("sc1a1drY"), "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, selected = sc1def$dimred[2])
          )
        )
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
        6, style = "border-right: 2px solid black", h4("Cell information"),
        fluidRow(
          column(
            6,
            selectInput(
              ns("sc1a1inp1"),
              "Cell information:",
              choices = sc1conf$UI,
              selected = sc1def$meta1
            ) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Cell information to colour cells by",
                content = c(
                  "Select cell information to colour cells",
                  "Categorical covariates have a fixed colour palette",
                  paste0(
                    "Continuous covariates are coloured in a ",
                    "Blue-Yellow-Red colour scheme, which can be ",
                    "changed in the plot controls"
                  )
                )
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
        
        fluidRow(column(12, uiOutput(ns("sc1a1oup1_ui")))),
        downloadButton(ns("sc1a1oup1_pdf"), "Download PDF"),
        downloadButton(ns("sc1a1oup1_png"), "Download PNG"),
        br(),
        
        div(
          style = "display:inline-block",
          numericInput(ns("sc1a1oup1_h"), "PDF / PNG height:", width = "138px", min = 4, max = 20, value = 6, step = 0.5)
        ),
        div(
          style = "display:inline-block",
          numericInput(ns("sc1a1oup1_w"), "PDF / PNG width:", width = "138px", min = 4, max = 20, value = 8, step = 0.5)
        ),
        br(),
        
        actionButton(ns("sc1a1tog9"), "Show cell numbers / statistics"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a1tog9")),
          h4("Cell numbers / statistics"),
          radioButtons(ns("sc1a1splt"), "Split continuous cell info into:",
                       choices = c("Quartile", "Decile"),
                       selected = "Decile", inline = TRUE),
          dataTableOutput(ns("sc1a1_dt"))
        )
      ),
      
      column(
        6, h4("Gene expression"),
        fluidRow(
          column(
            6,
            selectInput(ns("sc1a1inp2"), "Gene name:", choices = NULL) %>%
              helper(
                type = "inline", size = "m", fade = TRUE,
                title = "Gene expression to colour cells by",
                content = c(
                  "Select gene to colour cells by gene expression",
                  paste0(
                    "Gene expression are coloured in a ",
                    "White-Red colour scheme which can be ",
                    "changed in the plot controls"
                  )
                )
              )
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
        
        fluidRow(column(12, uiOutput(ns("sc1a1oup2_ui")))),
        downloadButton(ns("sc1a1oup2_pdf"), "Download PDF"),
        downloadButton(ns("sc1a1oup2_png"), "Download PNG"),
        br(),

        div(
          style = "display:inline-block",
          numericInput(ns("sc1a1oup2_h"), "PDF / PNG height:", width = "138px", min = 4, max = 20, value = 6, step = 0.5)
        ),
        div(
          style = "display:inline-block",
          numericInput(ns("sc1a1oup2_w"), "PDF / PNG width:", width = "138px", min = 4, max = 20, value = 8, step = 0.5)
        ),
        br(),
        
        ######################## ???MARKERGENES
        actionButton(ns("sc1a1tog10"), "Show Marker Genes Per Cluster"),
        
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a1tog10")),
          
          h4("Marker Genes"),
          
          radioButtons(
            ns("sc1a1splt_test"),
            "Order top genes by:",
            choices = c(
              "logFC and Adj Pval (Wilcox)",
              "AUC (ROC)",
              "Average Expression",
              "% Expression In Cluster"
            ),
            selected = "AUC (ROC)",
            inline = TRUE
          ),
          
          selectInput(
            ns("resolution"),
            "Clustering resolution:",
            choices = sc1conf$UI[grep("res", sc1conf$UI)],
            selected = sc1conf$UI[grep("res", sc1conf$UI)][1],
            multiple = FALSE
          ),
          
          checkboxInput(ns("show_all"), "Show all genes", value = FALSE),
          
          conditionalPanel(
            condition = sprintf("!input['%s']", ns("show_all")),
            sliderInput(ns("top"), "Number of genes per cluster", min = 1, max = 50, value = 10, step = 1)
          ),
          
          dataTableOutput(ns("sc1a1_dtmarkers"))
        )
        
        ##############################
      )
    )
  )
}

############################################### Server #################################################

scDRnum_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    ### For all tags and Server-side selectize 
    observe_helpers() 
    optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
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
                         choices = c(sc1conf[is.na(fID)]$UI,names(sc1gene)), 
                         selected = sc1conf[is.na(fID)]$UI[1], options = list( 
                           maxOptions = length(sc1conf[is.na(fID)]$UI) + 3, 
                           create = TRUE, persist = TRUE, render = I(optCrt))) 
    
    
    # Needed for plot sizes (ShinyCell defaults)
    # If pList exists globally (from sourced ShinyCell helpers), this does nothing.
    if (!exists("pList", inherits = TRUE)) {
      pList <<- c(Small = "350px", Medium = "550px", Large = "750px")
    }
    
    # Required: scDRcell and scDRgene must exist (from other sourced ShinyCell modules)
    if (!exists("scDRcell", inherits = TRUE)) {
      stop("scDRcell() is not available. Source the ShinyCell function file that defines scDRcell (often scDRcell.R).")
    }
    if (!exists("scDRgene", inherits = TRUE)) {
      stop("scDRgene() is not available. Source the ShinyCell function file that defines scDRgene (often scDRgene.R).")
    }
    
    output$sc1a1sub1_ui <- renderUI({
      req(input$sc1a1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(
        ns("sc1a1sub2"),
        "Select which cells to show",
        inline = TRUE,
        choices = sub,
        selected = sub
      )
    })
    
    observeEvent(input$sc1a1sub1non, {
      req(input$sc1a1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(
        session,
        inputId = "sc1a1sub2",
        label = "Select which cells to show",
        choices = sub,
        selected = NULL,
        inline = TRUE
      )
    })
    
    observeEvent(input$sc1a1sub1all, {
      req(input$sc1a1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(
        session,
        inputId = "sc1a1sub2",
        label = "Select which cells to show",
        choices = sub,
        selected = sub,
        inline = TRUE
      )
    })
    
    # Render plot 1 (cell info) and provide the UI container
    output$sc1a1oup1 <- renderPlot({
      req(input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1)
      scDRcell(
        sc1conf, sc1meta,
        input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,
        input$sc1a1sub1, input$sc1a1sub2,
        input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,
        input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1
      )
    })
    
    output$sc1a1oup1_ui <- renderUI({
      req(input$sc1a1psz)
      plotOutput(ns("sc1a1oup1"), height = pList[input$sc1a1psz])
    })
    
    output$sc1a1oup1_pdf <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1a1drX, "_", input$sc1a1drY, "_", input$sc1a1inp1, ".pdf")
      },
      content = function(file) {
        ggsave(
          file,
          device = "pdf",
          height = input$sc1a1oup1_h,
          width  = input$sc1a1oup1_w,
          useDingbats = FALSE,
          plot = scDRcell(
            sc1conf, sc1meta,
            input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,
            input$sc1a1sub1, input$sc1a1sub2,
            input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,
            input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1
          )
        )
      }
    )
    
    output$sc1a1oup1_png <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1a1drX, "_", input$sc1a1drY, "_", input$sc1a1inp1, ".png")
      },
      content = function(file) {
        ggsave(
          file,
          device = "png",
          height = input$sc1a1oup1_h,
          width  = input$sc1a1oup1_w,
          plot = scDRcell(
            sc1conf, sc1meta,
            input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,
            input$sc1a1sub1, input$sc1a1sub2,
            input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,
            input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1
          )
        )
      }
    )
    
    # Render plot 2 (gene expression) and provide the UI container
    output$sc1a1oup2 <- renderPlot({
      req(input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2)
      scDRgene(
        sc1conf, sc1meta,
        input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,
        input$sc1a1sub1, input$sc1a1sub2,
        file.path(dir_inputs, "sc1gexpr.h5"),
        sc1gene,
        input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2,
        input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt
      )
    })
    
    output$sc1a1oup2_ui <- renderUI({
      req(input$sc1a1psz)
      plotOutput(ns("sc1a1oup2"), height = pList[input$sc1a1psz])
    })
    
    output$sc1a1oup2_pdf <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1a1drX, "_", input$sc1a1drY, "_", input$sc1a1inp2, ".pdf")
      },
      content = function(file) {
        ggsave(
          file,
          device = "pdf",
          height = input$sc1a1oup2_h,
          width  = input$sc1a1oup2_w,
          useDingbats = FALSE,
          plot = scDRgene(
            sc1conf, sc1meta,
            input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,
            input$sc1a1sub1, input$sc1a1sub2,
            file.path(dir_inputs, "sc1gexpr.h5"),
            sc1gene,
            input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2,
            input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt
          )
        )
      }
    )
    
    output$sc1a1oup2_png <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1a1drX, "_", input$sc1a1drY, "_", input$sc1a1inp2, ".png")
      },
      content = function(file) {
        ggsave(
          file,
          device = "png",
          height = input$sc1a1oup2_h,
          width  = input$sc1a1oup2_w,
          plot = scDRgene(
            sc1conf, sc1meta,
            input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,
            input$sc1a1sub1, input$sc1a1sub2,
            file.path(dir_inputs, "sc1gexpr.h5"),
            sc1gene,
            input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2,
            input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt
          )
        )
      }
    )
    
    # Table (already present)
    output$sc1a1_dt <- renderDataTable({
      ggData <- scDRnum(
        sc1conf, sc1meta,
        input$sc1a1inp1, input$sc1a1inp2,
        input$sc1a1sub1, input$sc1a1sub2,
        file.path(dir_inputs, "sc1gexpr.h5"),
        sc1gene,
        input$sc1a1splt
      )
      
      datatable(
        ggData,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))
      ) %>%
        formatRound(columns = c("pctExpress"), digits = 2)
    })
    ################################ ???
    # Table (Marker Genes)
   
    output$sc1a1_dtmarkers <-  renderDT(server = FALSE, {
        
        req(markers_list)
     
        resolution_selection <- paste0(input$resolution)
        top_selection <- input$top
        
        ds <- arrow::open_dataset(markers_list)
        
        if (isTRUE(input$show_all)) {
          df <- ds |>
            dplyr::filter(annotation == resolution_selection) |>
            dplyr::collect()
          
          DT::datatable(
            df,
            extensions = c("Buttons"),
            options = list(
              dom = "Bfrtip",
              buttons = list(
                list(
                  extend = "csv",
                  text = "Download Full Results",
                  filename = paste0("MarkersList_", resolution_selection),
                  exportOptions = list(modifier = list(page = "all"))
                )
              )
            )
          )
        } else {
          top_gene <- as.integer(input$top)
          rank_by_selection<-input$sc1a1splt_test
          
          observeEvent(markers_list, {
            ds <- arrow::open_dataset(markers_list)
            message("marker columns: ", paste(names(ds), collapse = ", "))
          })
          
          if (rank_by_selection == "logFC and Adj Pval (Wilcox)") {
            
            df <- ds |>
              dplyr::filter(annotation == resolution_selection) |>
              dplyr::select(feature, group, logFC, padj) |>
              dplyr::collect()
            
            shiny::validate(need(all(c("padj", "logFC", "group") %in% names(df)),
                          "Expected columns not found. Check your parquet columns (padj, logFC, group)."))
            
            df <- df |>
              dplyr::filter(padj < 0.05) |>
              dplyr::group_by(group) |>
              dplyr::arrange(dplyr::desc(logFC), .by_group = TRUE) |>
              dplyr::slice_head(n = top_gene)
            
          } else if (rank_by_selection == "AUC (ROC)") {
            
            df <- ds |>
              dplyr::filter(annotation == resolution_selection) |>
              dplyr::select(feature, group, auc, pct_in, pct_out) |>
              dplyr::collect()
            
            shiny::validate(need(all(c("auc", "group") %in% names(df)),
                          "Expected columns not found. Check your parquet columns (auc, group)."))
            
            df <- df |>
              dplyr::group_by(group) |>
              dplyr::arrange(dplyr::desc(auc), .by_group = TRUE) |>
              dplyr::slice_head(n = top_gene)
            
          } else if (rank_by_selection == "Average Expression") {
            
            df <- ds |>
              dplyr::filter(annotation == resolution_selection) |>
              dplyr::select(feature, group, avgExpr, pct_in, pct_out) |>
              dplyr::collect()
            
            shiny::validate(need(all(c("avgExpr", "group") %in% names(df)),
                          "Expected columns not found. Check your parquet columns (avgExpr, group)."))
            
            df <- df |>
              dplyr::group_by(group) |>
              dplyr::arrange(dplyr::desc(avgExpr), .by_group = TRUE) |>
              dplyr::slice_head(n = top_gene)
            
          } else if (rank_by_selection == "% Expression In Cluster") {
            
            df <- ds |>
              dplyr::filter(annotation == resolution_selection) |>
              dplyr::select(feature, group, pct_in, pct_out) |>
              dplyr::collect()
            
            shiny::validate(need(all(c("pct_in", "group") %in% names(df)),
                          "Expected columns not found. Check your parquet columns (pct_in, group)."))
            
            df <- df |>
              dplyr::group_by(group) |>
              dplyr::arrange(dplyr::desc(pct_in), .by_group = TRUE) |>
              dplyr::slice_head(n = top_gene)
            
          } else {
            shiny::validate(need(FALSE, paste0("Unknown ranking option: ", rank_by_selection)))
          }
          
          shiny::validate(need(!is.null(df) && nrow(df) > 0, "No rows after filtering. Try a different resolution or relax filters."))
          
          DT::datatable(
            df,
            extensions = c("Buttons"),
            options = list(
              dom = "Bfrtip",
              buttons = list(
                list(
                  extend = "csv",
                  text = "Download Top Genes",
                  filename = paste0("SubsetMarkersList_top",top_gene,"_", resolution_selection),
                  exportOptions = list(modifier = list(page = "all"))
                )
              )
            )
          )
        }
      })
      
    ################################ 
  })}


############################################### Registration #################################################

register_tab(
  id     = "cellinfo_geneexpr",
  title  = "CellInfo vs GeneExpr",
  ui     = scDRnum_ui,
  server = scDRnum_server
)
