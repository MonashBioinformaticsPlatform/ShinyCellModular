### Functions for tab A1 
getGsc1a1inp1 <- reactive({ 
  req(gsub("^Assay: ", "", input$sc1a1ass1)) 
  if(gsub("^Assay: ", "", input$sc1a1ass1) == "Cell Information"){ 
    res <- sc1conf$UI 
    resDef <- sc1def$meta1 
    resLen <- length(res) 
  } else { 
    res <- names(sc1gene[[gsub("^Assay: ", "", input$sc1a1ass1)]]) 
    resDef <- sc1def$gene1[[gsub("^Assay: ", "", input$sc1a1ass1)]] 
    resLen <- 7 
  } 
  return(list(res,resDef,resLen))
})
observeEvent(gsub("^Assay: ", "", input$sc1a1ass1), {
  updateSelectizeInput(session, "sc1a1inp1", choices = getGsc1a1inp1()[[1]], 
                       server = TRUE, selected = getGsc1a1inp1()[[2]], options = list( 
                         maxOptions = getGsc1a1inp1()[[3]], create = TRUE, 
                         persist = TRUE, render = I(optCrt))) 
})
output$sc1a1sub1.ui <- renderUI({ 
  sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
  checkboxGroupInput("sc1a1sub2", "Select which cells to show", inline = TRUE, 
                     choices = sub, selected = sub) 
}) 
observeEvent(input$sc1a1sub1non, { 
  sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
  updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                           choices = sub, selected = NULL, inline = TRUE) 
}) 
observeEvent(input$sc1a1sub1all, { 
  sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
  updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                           choices = sub, selected = sub, inline = TRUE) 
}) 

sc1a1oup1xy <- reactiveValues(x = NULL, y = NULL)
observe({
  brush <- input$sc1a1inp1.br
  if (!is.null(brush)) {
    sc1a1oup1xy$x <- c(brush$xmin, brush$xmax)
    sc1a1oup1xy$y <- c(brush$ymin, brush$ymax)
  } else {
    sc1a1oup1xy$x <- NULL; sc1a1oup1xy$y <- NULL
  }
})
sc1a1oup1br <- reactive({
  sc2Ddimr(sc1conf, sc1meta, sc1dimr, input$sc1a1dr, input$sc1a1inp1,
           "sc1assay_", sc1gene, input$sc1a1ass1, input$sc1a1sub1, input$sc1a1sub2,  
           input$sc1a1min1, input$sc1a1max1, input$sc1a1siz/2, input$sc1a1ord1,
           cList[[input$sc1a1col1]], sList[input$sc1a1fsz]/2, 
           input$sc1a1asp, FALSE, FALSE) 
})
output$sc1a1oup1.br <- renderPlot({sc1a1oup1br() + theme(legend.position = "none")}) 

sc1a1oup1 <- reactive({
  sc2Ddimr(sc1conf, sc1meta, sc1dimr, input$sc1a1dr, input$sc1a1inp1,
           "sc1assay_", sc1gene, input$sc1a1ass1, input$sc1a1sub1, input$sc1a1sub2,  
           input$sc1a1min1, input$sc1a1max1, input$sc1a1siz, input$sc1a1ord1,
           cList[[input$sc1a1col1]], sList[input$sc1a1fsz], 
           input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1)
})
output$sc1a1oup1 <- renderPlot({
  if(is.null(sc1a1oup1xy$x[1])){
    sc1a1oup1() + theme(legend.position = "none")
  } else {
    sc1a1oup1() + theme(legend.position = "none") + 
      scale_x_continuous(limits = sc1a1oup1xy$x, expand = c(0, 0)) + 
      scale_y_continuous(limits = sc1a1oup1xy$y, expand = c(0, 0)) 
  }
})
output$sc1a1oup1.ui <- renderUI({plotOutput("sc1a1oup1", height = pList[input$sc1a1psz])})
output$sc1a1oup1.dl <- downloadHandler(
  filename = function() { paste0("sc1",input$sc1a1dr,"_",input$sc1a1inp1,".",input$sc1a1oup1.f) },
  content = function(file) { 
    if(is.null(sc1a1oup1xy$x[1])){
      ggsav(file, height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, 
            plot = sc1a1oup1() + theme(legend.position = "none"))
    } else {
      ggsav(file, height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, 
            plot = sc1a1oup1() + theme(legend.position = "none") + 
              scale_x_continuous(limits = sc1a1oup1xy$x, expand = c(0, 0)) + 
              scale_y_continuous(limits = sc1a1oup1xy$y, expand = c(0, 0)))
    }
  })
output$sc1a1oup3 <- renderPlot({grid.newpage(); grid.draw(g_legend(sc1a1oup1()))})
output$sc1a1oup3.ui <- renderUI({plotOutput("sc1a1oup3", height = 72*convertHeight(
  grobHeight(g_legend(sc1a1oup1())), unitTo="in", valueOnly=TRUE) + 50)})
output$sc1a1oup3.dl <- downloadHandler(
  filename = function() { paste0("sc1",input$sc1a1dr,"_",input$sc1a1inp1,"_leg.",input$sc1a1oup3.f) },
  content = function(file) { 
    grid.newpage(); grid.draw(g_legend(sc1a1oup1()))
    ggsav(file, height = input$sc1a1oup3.h, width = input$sc1a1oup3.w, plot = grid.grab())
  }) 

output$sc1a1.dt <- renderDataTable({
  ggData = sc2Dnum(sc1conf, sc1meta, sc1dimr, input$sc1a1dr, sc1a1oup1xy$x, sc1a1oup1xy$y, input$sc1a1inp1,
                   "sc1assay_", sc1gene, input$sc1a1ass1, input$sc1a1sub1, input$sc1a1sub2, input$sc1a1splt)
  datatable(ggData, rownames = FALSE, extensions = "Buttons",
            options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
    formatRound(columns = c("pctZoom"), digits = 2)
})
