############################################### Functions ############################################

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

################################################# UI #################################################
scDRgene_ui <- function(id, sc1conf, sc1def) {
  ns <- NS(id)
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc1a3drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
            selectInput("sc1a3drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                        selected = sc1def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc1a3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1a3togL % 2 == 1", 
          selectInput("sc1a3sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1a3sub1.ui"), 
          actionButton("sc1a3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1a3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc1a3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1a3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc1a3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc1a3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc1a3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc1a3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc1a3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a3tog1 % 2 == 1", 
              radioButtons("sc1a3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc1a3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a3oup1.ui"))), 
        downloadButton("sc1a3oup1.pdf", "Download PDF"), 
        downloadButton("sc1a3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a3tog2 % 2 == 1", 
              radioButtons("sc1a3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc1a3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a3oup2.ui"))), 
        downloadButton("sc1a3oup2.pdf", "Download PDF"), 
        downloadButton("sc1a3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 
}

############################################## Server ################################################

scDRgene_server <- function(id, sc1conf, sc1meta, sc1gene, dir_inputs) {
  moduleServer(id, function(input, output, session) {

        # you can use ns() for things you create inside renderUI
        ns <- session$ns

      output$sc1a1oup2 <- renderPlot({ 
        scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                input$sc1a1sub1, input$sc1a1sub2, 
                paste0(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) 
      }) 
      output$sc1a1oup2.ui <- renderUI({ 
        plotOutput("sc1a1oup2", height = pList[input$sc1a1psz]) 
      }) 
      output$sc1a1oup2.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                      input$sc1a1inp2,".pdf") }, 
        content = function(file) { ggsave( 
          file, device = "pdf", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, useDingbats = FALSE, 
          plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                          input$sc1a1sub1, input$sc1a1sub2, 
                          paste0(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                          input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                          input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
      }) 
      output$sc1a1oup2.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                      input$sc1a1inp2,".png") }, 
        content = function(file) { ggsave( 
          file, device = "png", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, 
          plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                          input$sc1a1sub1, input$sc1a1sub2, 
                          paste0(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                          input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                          input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
      }) 

      ### Plots for tab a3 
      output$sc1a3sub1.ui <- renderUI({ 
        sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
        checkboxGroupInput("sc1a3sub2", "Select which cells to show", inline = TRUE, 
                          choices = sub, selected = sub) 
      }) 
      observeEvent(input$sc1a3sub1non, { 
        sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                                choices = sub, selected = NULL, inline = TRUE) 
      }) 
      observeEvent(input$sc1a3sub1all, { 
        sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                                choices = sub, selected = sub, inline = TRUE) 
      }) 
      output$sc1a3oup1 <- renderPlot({ 
        scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                input$sc1a3sub1, input$sc1a3sub2, 
                paste0(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
      }) 
      output$sc1a3oup1.ui <- renderUI({ 
        plotOutput("sc1a3oup1", height = pList[input$sc1a3psz]) 
      }) 
      output$sc1a3oup1.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                      input$sc1a3inp1,".pdf") }, 
        content = function(file) { ggsave( 
          file, device = "pdf", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, useDingbats = FALSE, 
          plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                          input$sc1a3sub1, input$sc1a3sub2, 
                          paste0(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                          input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                          input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
      }) 
      output$sc1a3oup1.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                      input$sc1a3inp1,".png") }, 
        content = function(file) { ggsave( 
          file, device = "png", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, 
          plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                          input$sc1a3sub1, input$sc1a3sub2, 
                          paste0(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                          input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                          input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
      }) 
      
      output$sc1a3oup2 <- renderPlot({ 
        scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                input$sc1a3sub1, input$sc1a3sub2, 
                paste0(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
      }) 
      output$sc1a3oup2.ui <- renderUI({ 
        plotOutput("sc1a3oup2", height = pList[input$sc1a3psz]) 
      }) 
      output$sc1a3oup2.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                      input$sc1a3inp2,".pdf") }, 
        content = function(file) { ggsave( 
          file, device = "pdf", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, useDingbats = FALSE, 
          plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                          input$sc1a3sub1, input$sc1a3sub2, 
                          paste0(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                          input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                          input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
      }) 
      output$sc1a3oup2.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                      input$sc1a3inp2,".png") }, 
        content = function(file) { ggsave( 
          file, device = "png", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, 
          plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                          input$sc1a3sub1, input$sc1a3sub2, 
                          paste0(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                          input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                          input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
      }) 
  })  # closes moduleServer
} 
     

     