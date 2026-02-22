
ui <- navbarPage(
  "",
  tabPanel(
    title = "App",
    fluidPage( theme = shinytheme("cerulean"),
               tags$head(
                 tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}")),
                 tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))
               ),
               titlePanel(app_title),
               do.call(navbarPage, c(list(title = NULL), tab_panels)),
               tags$hr(),
               tags$p(
                 style = "font-size: 90%%; color: #666;",
                 em(
                   "This application was generated using ShinyCellPlus. Monash Genomics and Bioinformatics Platform",
                   "Tabs are dynamically loaded from modular components.ShinyCellPlus: ShinyCell Package Customized by MGBP v.1 Date: Jan 2026"
                 )
               ),
               br(), br(), br(), br(), br()
    )
  )
)

# This tab is base in the orginal "cellinfo_geneexpr" in ShinyCell 
# id     = "cellinfo_geneexpr",
# title  = "CellInfo vs GeneExpr",

############################################### Functions ############################################
actionButton(ns("sc1a1tog10"), "Show Marker Genes Per Cluster"),
conditionalPanel(
  condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a1tog10")),
  h4("Marker Genes"),
  radioButtons(ns("sc1a1splt_test"), "Order top genes by:",
               choices = c("logFC and AdjPvalue (Wilcox)", "AUC (ROC)","Average Expression","Percentage of Expression IN Cluster"),
               selected = "AUC", inline = TRUE),
  selectInput("resolution","Clustering resolution:",
              choices= sc1conf$UI[grep("res",sc1conf$UI)],
              selected  = sc1conf$UI[grep("res",sc1conf$UI)][1],multiple = FALSE),
  condition = sprintf("!input['%s']", ns("show_all")),
  sliderInput(ns("top"), "Number of genes per cluster", min = 1, max = 50, value = 10, step = 1),
  dataTableOutput(ns("sc1a1_dtmarkers"))
)

scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                    inpH5, inpGene, inpsplt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("group", "sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Split inp1 if necessary 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    if(inpsplt == "Quartile"){nBk = 4} 
    if(inpsplt == "Decile"){nBk = 10} 
    ggData$group = cut(ggData$group, breaks = nBk) 
  } 
  
  # Actual data.table 
  ggData$express = FALSE 
  ggData[val2 > 0]$express = TRUE 
  ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
  ggData = ggData[, .(nCells = .N), by = "group"] 
  ggData = ggData1[ggData, on = "group"] 
  ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
  ggData[is.na(nExpress)]$nExpress = 0 
  ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
  ggData = ggData[order(group)] 
  colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
  return(ggData) 
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
        actionButton(ns("sc1a1togL"), "Toggle to subset cells"),
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
        actionButton(ns("sc1a1tog0"), "Toggle graphics controls"),
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
            actionButton(ns("sc1a1tog1"), "Toggle plot controls"),
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
        
        actionButton(ns("sc1a1tog9"), "Toggle to show cell numbers / statistics"),
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
            actionButton(ns("sc1a1tog2"), "Toggle plot controls"),
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
        )
      )
    )
  )
}

############################################### Server #################################################


scDRnum_server <- function(id, sc1conf, sc1meta, sc1gene, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    output$sc1a1sub1_ui <- renderUI({
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
  })
}

############################################### Registration #################################################

register_tab(
  id     = "cellinfo_geneexpr",
  title  = "CellInfo vs GeneExpr",
  ui     = scDRnum_ui,
  server = scDRnum_server
)





# This tab is base in the orginal "cellinfo_geneexpr" in ShinyCell 
# id     = "cellinfo_geneexpr",
# title  = "CellInfo vs GeneExpr",

############################################### Functions ############################################


scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                    inpH5, inpGene, inpsplt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("group", "sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Split inp1 if necessary 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    if(inpsplt == "Quartile"){nBk = 4} 
    if(inpsplt == "Decile"){nBk = 10} 
    ggData$group = cut(ggData$group, breaks = nBk) 
  } 
  
  # Actual data.table 
  ggData$express = FALSE 
  ggData[val2 > 0]$express = TRUE 
  ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
  ggData = ggData[, .(nCells = .N), by = "group"] 
  ggData = ggData1[ggData, on = "group"] 
  ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
  ggData[is.na(nExpress)]$nExpress = 0 
  ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
  ggData = ggData[order(group)] 
  colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
  return(ggData) 
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
        actionButton(ns("sc1a1togL"), "Toggle to subset cells"),
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
        actionButton(ns("sc1a1tog0"), "Toggle graphics controls"),
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
            actionButton(ns("sc1a1tog1"), "Toggle plot controls"),
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
        
        actionButton(ns("sc1a1tog9"), "Toggle to show cell numbers / statistics"),
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
            actionButton(ns("sc1a1tog2"), "Toggle plot controls"),
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
        )
      )
    )
  )
}

############################################### Server #################################################


scDRnum_server <- function(id, sc1conf, sc1meta, sc1gene, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    output$sc1a1sub1_ui <- renderUI({
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
  })
}

############################################### Registration #################################################

register_tab(
  id     = "cellinfo_geneexpr",
  title  = "CellInfo vs GeneExpr",
  ui     = scDRnum_ui,
  server = scDRnum_server
)