# This tab is base in the orginal "genecoex" in ShinyCell
# id     = "genecoex",
# title  = "Gene Coexpression",

##################################################### Functions #########################################################
# Plot gene coexpression on dimred 
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){ 
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22 
  oup = oup / (xy*xy) 
  return(oup) 
} 

scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Map colours 
  ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1)) 
  ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2)) 
  ggData$v0 = ggData$v1 + ggData$v2 
  ggData = gg[ggData, on = c("v1", "v2")] 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(v0)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-v0)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16, color = ggData$cMix) + 
    xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) + 
    scale_color_gradientn(inp1, colours = cList[[1]]) + 
    guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){ 
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Actual ggplot 
  ggOut = ggplot(gg, aes(v1, v2)) + 
    geom_tile(fill = gg$cMix) + 
    xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) + 
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    sctheme(base_size = sList[inpfsz], XYval = TRUE) 
  return(ggOut) 
} 
 
scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE] 
  colnames(ggData) = c("sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Actual data.table 
  ggData$express = "none" 
  ggData[val1 > 0]$express = inp1 
  ggData[val2 > 0]$express = inp2 
  ggData[val1 > 0 & val2 > 0]$express = "both" 
  ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none"))) 
  ggData = ggData[, .(nCells = .N), by = "express"] 
  ggData$percent = 100 * ggData$nCells / sum(ggData$nCells) 
  ggData = ggData[order(express)] 
  colnames(ggData)[1] = "expression > 0" 
  return(ggData) 
} 

##################################################### UI #########################################################

scDRcoex_ui <- function(id, sc1conf, sc1def){
  ns <- NS(id)
  
  tabPanel( 
    HTML("Gene coexpression"), 
    h4("Coexpression of two genes on reduced dimensions"), 
    "In this tab, users can visualise the coexpression of two genes ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput(ns("sc1b2drX"), "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                            selected = sc1def$dimred[1]), 
            selectInput(ns("sc1b2drY"), "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                        selected = sc1def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton(ns("sc1b2togL"), "Filter Cells"), 
        conditionalPanel( 
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1b2togL")),
          selectInput(ns("sc1b2sub1"), "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput(ns("sc1b2sub1.ui")), 
          actionButton(ns("sc1b2sub1all"), "Select all groups", class = "btn btn-primary"), 
          actionButton(ns("sc1b2sub1non"), "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton(ns("sc1b2tog0"), "Customize plot"), 
        conditionalPanel( 
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1b2tog0")),
          fluidRow( 
            column( 
              6, sliderInput(ns("sc1b2siz"), "Point size:", 
                              min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons(ns("sc1b2psz"), "Plot size:", 
                            choices = c("Small", "Medium", "Large"), 
                            selected = "Medium", inline = TRUE), 
              radioButtons(ns("sc1b2fsz"), "Font size:", 
                            choices = c("Small", "Medium", "Large"), 
                            selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons(ns("sc1b2asp"), "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput(ns("sc1b2txt"), "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", h4("Gene Expression"), 
        selectInput(ns("sc1b2inp1"), "Gene 1:", choices=NULL) %>%  
          helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                  "White-Red colour scheme which can be ", 
                                  "changed in the plot controls"))), 
        selectInput(ns("sc1b2inp2"), "Gene 2:", choices=NULL) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                  title = "Gene expression to colour cells by", 
                  content = c("Select gene to colour cells by gene expression", 
                              paste0("- Gene expression are coloured in a ", 
                                    "White-Blue colour scheme which can be ", 
                                    "changed in the plot controls"))), 
        actionButton(ns("sc1b2tog1"), "Customize plot"), 
        conditionalPanel( 
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1b2tog1")),
          radioButtons(ns("sc1b2col1"), "Colour:", 
                        choices = c("Red (Gene1); Blue (Gene2)", 
                                    "Orange (Gene1); Blue (Gene2)", 
                                    "Red (Gene1); Green (Gene2)", 
                                    "Green (Gene1); Blue (Gene2)"), 
                        selected = "Red (Gene1); Blue (Gene2)"), 
          radioButtons(ns("sc1b2ord1"), "Plot order:", 
                        choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                        selected = "Max-1st", inline = TRUE) 
        ) 
      ), # End of column (6 space) 
      column( 
        6, style="border-right: 2px solid black", 
        uiOutput(ns("sc1b2oup1.ui")), 
        downloadButton(ns("sc1b2oup1.pdf"), "Download PDF"), 
        downloadButton(ns("sc1b2oup1.png"), "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput(ns("sc1b2oup1.h"), "PDF / PNG height:", width = "138px", 
                          min = 4, max = 20, value = 8, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1b2oup1.w", "PDF / PNG width:", width = "138px", 
                          min = 4, max = 20, value = 10, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        3, uiOutput(ns("sc1b2oup2.ui")), 
        downloadButton(ns("sc1b2oup2.pdf"), "Download PDF"), 
        downloadButton(ns("sc1b2oup2.png"), "Download PNG"), 
        br(), h4("Cell numbers"), 
        dataTableOutput(ns("sc1b2.dt"))
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )
}


##################################################### Server #########################################################

scDRcoex_server <-  function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
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
      
      
  
  ### Plots for tab b2 
  output$sc1b2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1b2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1b2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1b2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
     
  output$sc1b2oup1 <- renderPlot({ 
    scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,   
             input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
             file.path(dir_inputs,"sc1gexpr.h5"), sc1gene, 
             input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
             input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) 
  }) 
  output$sc1b2oup1.ui <- renderUI({ 
    plotOutput(ns("sc1b2oup1"), height = pList2[input$sc1b2psz]) 
  }) 
  output$sc1b2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp1,"_",input$sc1b2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                      input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                      file.path(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
  }) 
  output$sc1b2oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp1,"_",input$sc1b2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                      input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                      file.path(dir_inputs,"sc1gexpr.h5"), sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
  }) 
  output$sc1b2oup2 <- renderPlot({ 
    scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) 
  }) 
  output$sc1b2oup2.ui <- renderUI({ 
    plotOutput(ns("sc1b2oup2"), height = "300px") 
  }) 
  output$sc1b2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
  }) 
  output$sc1b2oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
  }) 
  output$sc1b2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(sc1conf, sc1meta, input$sc1b2inp1, input$sc1b2inp2, 
                         input$sc1b2sub1, input$sc1b2sub2, file.path(dir_inputs,"sc1gexpr.h5"), sc1gene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 

    })
}

############################################### Registration #################################################

register_tab(
  id     = "genecoex_multi",
  title  = "Gene Coexpression",
  ui     = scDRcoex_ui,
  server = scDRcoex_server
)


  