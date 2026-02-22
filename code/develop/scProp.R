

############################################### Functions ############################################
# Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                   inptyp, inpflp, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "grp", "sub") 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  ggData = ggData[, .(nCells = .N), by = c("X", "grp")] 
  ggData = ggData[, {tot = sum(nCells) 
                      .SD[,.(pctCells = 100 * sum(nCells) / tot, 
                             nCells = nCells), by = "grp"]}, by = "X"] 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$grp) 
  ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)] 
  ggData$grp = factor(ggData$grp, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "Proportion"){ 
    ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) + 
      geom_col() + ylab("Cell Proportion (%)") 
  } else { 
    ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) + 
      geom_col() + ylab("Number of Cells") 
  } 
  if(inpflp){ 
    ggOut = ggOut + coord_flip() 
  } 
  ggOut = ggOut + xlab(inp1) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) + 
    theme(legend.position = "right") 
  return(ggOut) 
} 
 
############################################## UI ################################################

scProp_ui <- function(id, sc1conf, sc1def) {
 
  tabPanel( 
    HTML("Proportion plot"), 
    h4("Proportion / cell numbers across different cell information"), 
    "In this tab, users can visualise the composition of single cells based on one discrete ", 
    "cell information across another discrete cell information. ",  
    "Usage examples include the library or cellcycle composition across clusters.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        selectInput("sc1c2inp1", "Cell information to plot (X-axis):", 
                    choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1def$grp2) %>%  
          helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to plot cells by",  
                content = c("Select categorical cell information to plot cells by", 
                            "- Plotted as the X-axis of the proportion plot")), 
        selectInput("sc1c2inp2", "Cell information to group / colour by:", 
                    choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1def$grp1) %>%  
          helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group / colour cells by", 
                content = c("Select categorical cell information to group / colour cells by", 
                            "- Proportion / cell numbers are shown in different colours")), 
        radioButtons("sc1c2typ", "Plot value:", 
                    choices = c("Proportion", "CellNumbers"), 
                    selected = "Proportion", inline = TRUE), 
        checkboxInput("sc1c2flp", "Flip X/Y", value = FALSE), 
        actionButton("sc1c2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1c2togL % 2 == 1", 
          selectInput("sc1c2sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1c2sub1.ui"), 
          actionButton("sc1c2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1c2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("sc1c2tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1c2tog % 2 == 1", 
          radioButtons("sc1c2psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
          radioButtons("sc1c2fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, uiOutput("sc1c2oup.ui"),  
            downloadButton("sc1c2oup.pdf", "Download PDF"),  
            downloadButton("sc1c2oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("sc1c2oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("sc1c2oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space)
}

############################################## Server ################################################
 
 
 scProp_server <- function(id, sc1conf, sc1def) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    ### Plots for tab c2 

    # UI elements
      output$sc1c2sub1.ui <- renderUI({ 
        sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
        checkboxGroupInput("sc1c2sub2", "Select which cells to show", inline = TRUE, 
                          choices = sub, selected = sub) 
      }) 
      observeEvent(input$sc1c2sub1non, { 
        sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                                choices = sub, selected = NULL, inline = TRUE) 
      }) 
      observeEvent(input$sc1c2sub1all, { 
        sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                                choices = sub, selected = sub, inline = TRUE) 
      }) 

    # Plot output

    output$sc1c2oup <- renderPlot({ 
      scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
            input$sc1c2sub1, input$sc1c2sub2, 
            input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) 
    }) 

    output$sc1c2oup.ui <- renderUI({ 
      plotOutput("sc1c2oup", height = pList2[input$sc1c2psz]) 
    }) 

    # Download handlers
    output$sc1c2oup.pdf <- downloadHandler( 
      filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                    input$sc1c2inp2,".pdf") }, 
      content = function(file) { ggsave( 
        file, device = "pdf", height = input$sc1c2oup.h, width = input$sc1c2oup.w, useDingbats = FALSE, 
        plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                      input$sc1c2sub1, input$sc1c2sub2, 
                      input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
      }) 
    output$sc1c2oup.png <- downloadHandler( 
      filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                    input$sc1c2inp2,".png") }, 
      content = function(file) { ggsave( 
        file, device = "png", height = input$sc1c2oup.h, width = input$sc1c2oup.w, 
        plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                      input$sc1c2sub1, input$sc1c2sub2, 
                      input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
      }) 
        
  })
} 