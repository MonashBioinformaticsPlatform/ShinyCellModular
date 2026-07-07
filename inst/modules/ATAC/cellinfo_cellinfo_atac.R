# CellInfo vs CellInfo (ATAC)
# This tab mirrors "cellinfo_cellinfo" (RNA), same logic, ATAC assay
# id     = "cellinfo_cellinfo_atac",
# title  = "CellInfo vs CellInfo (ATAC)",
########################################## Function ########################################################

# Plot cell information on dimred 
scDRcell_atac <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
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

########################################## UI ########################################################
cellinfo_cellinfo_atac_ui <- function(id, sc1conf_atac, sc1def_atac) {
  ns <- NS(id)
  
  tabPanel( 
    HTML("CellInfo vs CellInfo (ATAC)"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput(ns("sc1a2drX"), "X-axis:", choices = sc1conf_atac[dimred == TRUE]$UI, 
                            selected = sc1def_atac$dimred[1]), 
            selectInput(ns("sc1a2drY"), "Y-axis:", choices = sc1conf_atac[dimred == TRUE]$UI, 
                        selected = sc1def_atac$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton(ns("sc1a2togL"), "Filter cells"), 
        conditionalPanel( 
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2togL")),
          selectInput(ns("sc1a2sub1"), "Cell information to subset:", 
                      choices = sc1conf_atac[grp == TRUE]$UI, 
                      selected = sc1def_atac$grp1), 
          uiOutput(ns("sc1a2sub1.ui")), 
          actionButton(ns("sc1a2sub1all"), "Select all groups", class = "btn btn-primary"), 
          actionButton(ns("sc1a2sub1non"), "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton(ns("sc1a2tog0"), "Customize Aesthetics for Both Plots"), 
        conditionalPanel( 
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2tog0")) ,
          fluidRow( 
            column( 
              6, sliderInput(ns("sc1a2siz"), "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons(ns("sc1a2psz"), "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons(ns("sc1a2fsz"), "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons(ns("sc1a2asp"), "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput(ns("sc1a2txt"), "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput(ns("sc1a2inp1"), "Cell information:", 
                           choices = sc1conf_atac$UI, 
                           selected = sc1def_atac$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton(ns("sc1a2tog1"), "Customize plot"), 
            conditionalPanel( 
              condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2tog1")) ,
              radioButtons(ns("sc1a2col1"), "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons(ns("sc1a2ord1"), "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput(ns("sc1a2lab1"), "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput(ns("sc1a2oup1.ui")))), 
        downloadButton(ns("sc1a2oup1.pdf"), "Download PDF"), 
        downloadButton(ns("sc1a2oup1.png"), "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput(ns("sc1a2oup1.h"), "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput(ns("sc1a2oup1.w"), "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput(ns("sc1a2inp2"), "Cell information:", 
                           choices = sc1conf_atac$UI, 
                           selected = sc1def_atac$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton(ns("sc1a2tog2"), "Customize plot"), 
            conditionalPanel( 
              condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2tog2")), 
              radioButtons(ns("sc1a2col2"), "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons(ns("sc1a2ord2"), "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput(ns("sc1a2lab2"), "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput(ns("sc1a2oup2.ui"))), 
                 downloadButton(ns("sc1a2oup2.pdf"), "Download PDF"), 
                 downloadButton(ns("sc1a2oup2.png"), "Download PNG"), br(), 
                 div(style="display:inline-block", 
                     numericInput(ns("sc1a2oup2.h"), "PDF / PNG height:", width = "138px", 
                                  min = 4, max = 20, value = 6, step = 0.5)), 
                 div(style="display:inline-block", 
                     numericInput(ns("sc1a2oup2.w"), "PDF / PNG width:", width = "138px", 
                                  min = 4, max = 20, value = 8, step = 0.5)) 
        )
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )
}
########################################## server ########################################################

cellinfo_cellinfo_atac_server <- function(id, sc1conf_atac, sc1meta_atac, sc1def_atac, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    ### For all tags and Server-side selectize 
    observe_helpers() 
    
    ### Plots for tab a2 
    output$sc1a2sub1.ui <- renderUI({ 
      sub = strsplit(sc1conf_atac[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
      checkboxGroupInput(ns("sc1a2sub2"), "Select which cells to show", inline = TRUE, 
                         choices = sub, selected = sub) 
    }) 
    observeEvent(input$sc1a2sub1non, { 
      sub = strsplit(sc1conf_atac[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
      updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                               choices = sub, selected = NULL, inline = TRUE) 
    }) 
    observeEvent(input$sc1a2sub1all, { 
      sub = strsplit(sc1conf_atac[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
      updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                               choices = sub, selected = sub, inline = TRUE) 
    }) 
    output$sc1a2oup1 <- renderPlot({ 
      scDRcell_atac(sc1conf_atac, sc1meta_atac, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,  
               input$sc1a2sub1, input$sc1a2sub2, 
               input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1, 
               input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) 
    }) 
    output$sc1a2oup1.ui <- renderUI({ 
      plotOutput(ns("sc1a2oup1"), height = pList[input$sc1a2psz]) 
    }) 
    output$sc1a2oup1.pdf <- downloadHandler( 
      filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                     input$sc1a2inp1,".pdf") }, 
      content = function(file) { ggsave( 
        file, device = "pdf", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, useDingbats = FALSE, 
        plot = scDRcell_atac(sc1conf_atac, sc1meta_atac, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                        input$sc1a2sub1, input$sc1a2sub2, 
                        input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                        input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
      }) 
    output$sc1a2oup1.png <- downloadHandler( 
      filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                     input$sc1a2inp1,".png") }, 
      content = function(file) { ggsave( 
        file, device = "png", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, 
        plot = scDRcell_atac(sc1conf_atac, sc1meta_atac, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                        input$sc1a2sub1, input$sc1a2sub2, 
                        input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                        input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
      }) 
    
    output$sc1a2oup2 <- renderPlot({ 
      scDRcell_atac(sc1conf_atac, sc1meta_atac, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,  
               input$sc1a2sub1, input$sc1a2sub2, 
               input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2, 
               input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) 
    }) 
    output$sc1a2oup2.ui <- renderUI({ 
      plotOutput(ns("sc1a2oup2"), height = pList[input$sc1a2psz]) 
    }) 
    output$sc1a2oup2.pdf <- downloadHandler( 
      filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                     input$sc1a2inp2,".pdf") }, 
      content = function(file) { ggsave( 
        file, device = "pdf", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, useDingbats = FALSE, 
        plot = scDRcell_atac(sc1conf_atac, sc1meta_atac, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                        input$sc1a2sub1, input$sc1a2sub2, 
                        input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                        input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
      }) 
    output$sc1a2oup2.png <- downloadHandler( 
      filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                     input$sc1a2inp2,".png") }, 
      content = function(file) { ggsave( 
        file, device = "png", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, 
        plot = scDRcell_atac(sc1conf_atac, sc1meta_atac, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                        input$sc1a2sub1, input$sc1a2sub2, 
                        input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                        input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
      }) 
  }
  )
}



############################################### Registration #################################################
register_tab(
  id          = "cellinfo_cellinfo_atac",
  title       = "CellInfo vs CellInfo (ATAC)",
  ui          = cellinfo_cellinfo_atac_ui,
  server      = cellinfo_cellinfo_atac_server,
  author      = "Laura Perlaza-Jimenez",
  description = "Cell metadata visualisation against cell metadata, ATAC assay version of cellinfo_cellinfo",
  version     = "1.0",
  date        = "Jul 2026",
  source      = "Original ShinyCell, MGBP custom, ATAC adaptation",
  contact     = "laura.perlaza-jimenez@monash.edu"
)
