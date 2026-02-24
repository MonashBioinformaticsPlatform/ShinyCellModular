############################################### Functions ###########################################

# Get gene list
scGeneList <- function(inp, inpGene) {
  
  if (is.null(inp) || is.na(inp) || !nzchar(inp)) {
    return(data.table::data.table(gene = character(0), present = logical(0)))
  }
  
  # split on comma, semicolon, or newline (one or more)
  toks <- unlist(strsplit(inp, "[,;\\n\\r]+", perl = TRUE), use.names = FALSE)
  toks <- trimws(toks)
  toks <- toks[nzchar(toks)]
  
  geneList <- data.table::data.table(
    gene = unique(toks),
    present = TRUE
  )
  
  geneList[!gene %in% names(inpGene), present := FALSE]
  geneList[]
}
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt,
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol,
                       inpcols, inpfsz,
                       inpEngine = c("Classic ggplot", "FlexDotPlot"),
                       save = FALSE){
  
  inpEngine <- match.arg(inpEngine)
  
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  
  geneList = scGeneList(inp, inpGene)
  geneList = geneList[present == TRUE]
  shiny::validate(shiny::need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(shiny::need(nrow(geneList) > 1, "Please input at least 2 genes to plot!"))
  
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData = data.table()
  for(iGene in geneList$gene){
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) = c("sampleID", "sub")
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName = iGene
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=)))
    ggData = rbindlist(list(ggData, tmp))
  }
  h5file$close_all()
  
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% inpsub2]
  }
  shiny::validate(shiny::need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!"))
  
  ggData$val = expm1(ggData$val)
  ggData = ggData[, .(val = mean(val), prop = sum(val > 0) / length(sampleID)),
                  by = c("geneName", "grpBy")]
  ggData$val = log1p(ggData$val)
  
  colRange = range(ggData$val)
  if(inpScl){
    ggData[, val := scale(val), keyby = "geneName"]
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))
  }
  
  ggMat = dcast.data.table(ggData, geneName ~ grpBy, value.var = "val")
  tmp = ggMat$geneName
  ggMat = as.matrix(ggMat[, -1])
  rownames(ggMat) = tmp
  
  if(inpRow){
    hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat))))
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label)
  } else {
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene))
  }
  
  if(inpCol){
    hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label)
  }
  
  if(inpPlt == "Bubbleplot"){
    
    if(inpEngine == "FlexDotPlot"){
      
      shiny::validate(
        shiny::need(requireNamespace("FlexDotPlot", quietly = TRUE),
                    "FlexDotPlot engine selected but the FlexDotPlot package is not installed.")
      )
      
      dot_df <- as.data.frame(ggData[, .(grpBy, geneName, val, prop)])
      dot_df$grpBy   <- as.factor(dot_df$grpBy)
      dot_df$geneName <- as.factor(dot_df$geneName)
      
      res <- FlexDotPlot::dot_plot(
        data.to.plot = dot_df[, c("grpBy", "geneName")],
        size         = dot_df$prop,
        col          = dot_df$val,
        scale.by     = "radius",
        cols.use     = cList[[inpcols]],
        plot.legend  = TRUE,
        do.plot      = FALSE,
        do.return    = TRUE
      )
      
      return(res$plot)
    }
    
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) +
      geom_point() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_size_continuous("proportion", range = c(0, 8),
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(color = guide_colorbar(barwidth = 15)) +
      theme(axis.title = element_blank(), legend.box = "vertical")
    
  } else {
    
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) +
      geom_tile() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(fill = guide_colorbar(barwidth = 15)) +
      theme(axis.title = element_blank())
  }
  
  ggLeg = g_legend(ggOut)
  ggOut = ggOut + theme(legend.position = "none")
  
  if(!save){
    ggOut = grid.arrange(ggOut, ggLeg, heights = c(7,2), layout_matrix = rbind(c(1),c(2)))
  } else {
    ggOut = arrangeGrob(ggOut, ggLeg, heights = c(7,2), layout_matrix = rbind(c(1),c(2)))
  }
  
  return(ggOut)
}

################################################# UI #################################################

scBubbHeat_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
  tabPanel(
    HTML("Bubbleplot / Heatmap"),
    h4("Gene expression bubbleplot / heatmap"),
    "In this tab, users can visualise the gene expression patterns of ",
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(),
    "The normalised expression are averaged, log-transformed and then plotted.",
    br(), br(),
    
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        
        textAreaInput(
          ns("sc1d1inp"),
          HTML("List of gene names <br />
               (Max 50 genes, separated <br />
               by , or ; or newline):"),
          height = "200px",
          value = paste0(sc1def$genes, collapse = ", ")
        ),
        
        selectInput(
          ns("sc1d1grp"), "Group by:",
          choices = sc1conf[grp == TRUE]$UI,
          selected = sc1conf[grp == TRUE]$UI[1]
        ),
        
        radioButtons(
          ns("sc1d1plt"), "Plot type:",
          choices = c("Bubbleplot", "Heatmap"),
          selected = "Bubbleplot", inline = TRUE
        ),
        
        checkboxInput(ns("sc1d1scl"), "Scale gene expression", value = TRUE),
        checkboxInput(ns("sc1d1row"), "Cluster rows (genes)", value = TRUE),
        checkboxInput(ns("sc1d1col"), "Cluster columns (samples)", value = FALSE),
        
        br(),
        actionButton(ns("sc1d1togL"), "Filter Cells"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1d1togL")),
          selectInput(
            ns("sc1d1sub1"), "Cell information to subset:",
            choices = sc1conf[grp == TRUE]$UI,
            selected = sc1def$grp1
          ),
          uiOutput(ns("sc1d1sub1.ui")),
          actionButton(ns("sc1d1sub1all"), "Select all groups", class = "btn btn-primary"),
          actionButton(ns("sc1d1sub1non"), "Deselect all groups", class = "btn btn-primary")
        ),
        
        br(), br(),
        actionButton(ns("sc1d1tog"), "Customize plot"),
        
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1d1tog")),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Bubbleplot'", ns("sc1d1plt")),
            radioButtons(
              ns("sc1d1engine"),
              "Bubbleplot engine:",
              choices  = c("Classic ggplot", "FlexDotPlot"),
              selected = "Classic ggplot",
              inline   = TRUE
            )
          ),
          radioButtons(
            ns("sc1d1cols"), "Colour scheme:",
            choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
            selected = "Blue-Yellow-Red"
          ),
          radioButtons(
            ns("sc1d1psz"), "Plot size:",
            choices = c("Small", "Medium", "Large"),
            selected = "Medium", inline = TRUE
          ),
          radioButtons(
            ns("sc1d1fsz"), "Font size:",
            choices = c("Small", "Medium", "Large"),
            selected = "Medium", inline = TRUE
          )
        )
      ),
      
      column(
        9,
        h4(htmlOutput(ns("sc1d1oupTxt"))),
        uiOutput(ns("sc1d1oup.ui")),
        downloadButton(ns("sc1d1oup.pdf"), "Download PDF"),
        downloadButton(ns("sc1d1oup.png"), "Download PNG"),
        br(),
        div(style = "display:inline-block",
            numericInput(ns("sc1d1oup.h"), "PDF / PNG height:", width = "138px",
                         min = 4, max = 20, value = 10, step = 0.5)),
        div(style = "display:inline-block",
            numericInput(ns("sc1d1oup.w"), "PDF / PNG width:", width = "138px",
                         min = 4, max = 20, value = 10, step = 0.5))
      )
    )
  )
}

############################################## Server ################################################

scBubbHeat_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    observe_helpers()
    optCrt <- "{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
    
    updateSelectizeInput(session, "sc1a1inp2", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene1, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1a3inp1", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene1, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1a3inp2", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene2, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1b2inp1", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene1, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1b2inp2", choices = names(sc1gene), server = TRUE,
                         selected = sc1def$gene2, options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
    updateSelectizeInput(session, "sc1c1inp2", server = TRUE,
                         choices = c(sc1conf[is.na(fID)]$UI, names(sc1gene)),
                         selected = sc1conf[is.na(fID)]$UI[1], options = list(
                           maxOptions = length(sc1conf[is.na(fID)]$UI) + 3,
                           create = TRUE, persist = TRUE, render = I(optCrt)
                         ))
    
    if (!exists("pList3", inherits = TRUE)) {
      pList3 <<- c(Small = "450px", Medium = "650px", Large = "850px")
    }
    
    output$sc1d1sub1.ui <- renderUI({
      req(input$sc1d1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(ns("sc1d1sub2"), "Select which cells to show", inline = TRUE,
                         choices = sub, selected = sub)
    })
    
    observeEvent(input$sc1d1sub1non, {
      req(input$sc1d1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1d1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = NULL, inline = TRUE)
    })
    
    observeEvent(input$sc1d1sub1all, {
      req(input$sc1d1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1d1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = sub, inline = TRUE)
    })
    
    output$sc1d1oupTxt <- renderUI({
      geneList <- scGeneList(input$sc1d1inp, sc1gene)
      
      if (nrow(geneList) > 50) {
        return(HTML("More than 50 input genes. Please reduce the gene list."))
      }
      
      ok <- geneList[present == TRUE]
      notok <- geneList[present == FALSE]
      
      oup <- paste0(nrow(ok), " genes OK and will be plotted")
      if (nrow(notok) > 0) {
        oup <- paste0(
          oup, "<br/>",
          nrow(notok), " genes not found (",
          paste0(notok$gene, collapse = ", "),
          ")"
        )
      }
      HTML(oup)
    })
    
    
    output$sc1d1oup <- renderPlot({
      p <- scBubbHeat(
        sc1conf, sc1meta,
        input$sc1d1inp, input$sc1d1grp, input$sc1d1plt,
        input$sc1d1sub1, input$sc1d1sub2,
        paste0(dir_inputs,"sc1gexpr.h5"),
        sc1gene,
        input$sc1d1scl, input$sc1d1row, input$sc1d1col,
        input$sc1d1cols, input$sc1d1fsz,
        inpEngine = input$sc1d1engine
      )
      
      if (inherits(p, "grob") || inherits(p, "gtable")) {
        grid::grid.newpage()
        grid::grid.draw(p)
      } else {
        print(p)
      }
    })
    
    output$sc1d1oup.ui <- renderUI({
      plotOutput(ns("sc1d1oup"), height = pList3[input$sc1d1psz])
    })
    
    output$sc1d1oup.pdf <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1d1plt, "_", input$sc1d1grp, ".pdf")
      },
      content = function(file) {
        ggsave(
          file, device = "pdf",
          height = input$sc1d1oup.h, width = input$sc1d1oup.w,
          plot = scBubbHeat(
            sc1conf, sc1meta,
            input$sc1d1inp, input$sc1d1grp, input$sc1d1plt,
            input$sc1d1sub1, input$sc1d1sub2,
            paste0(dir_inputs,"sc1gexpr.h5"),
            sc1gene,
            input$sc1d1scl, input$sc1d1row, input$sc1d1col,
            input$sc1d1cols, input$sc1d1fsz,
            inpEngine = input$sc1d1engine,
            save = TRUE
          )
          )
      }
    )
    
    output$sc1d1oup.png <- downloadHandler(
      filename = function() {
        paste0("sc1", input$sc1d1plt, "_", input$sc1d1grp, ".png")
      },
      content = function(file) {
        ggsave(
          file, device = "png",
          height = input$sc1d1oup.h, width = input$sc1d1oup.w,
          plot = scBubbHeat(
            sc1conf, sc1meta,
            input$sc1d1inp,
            input$sc1d1grp,
            input$sc1d1plt,
            input$sc1d1sub1, input$sc1d1sub2,
            file.path(dir_inputs, "sc1gexpr.h5"),
            sc1gene,
            input$sc1d1scl, input$sc1d1row, input$sc1d1col,
            input$sc1d1cols, input$sc1d1fsz,
            save = TRUE
          )
        )
      }
    )
    
  })
}

############################################### Registration #################################################

register_tab(
  id     = "bubble_heatmap",
  title  = "Bubble Plot/ Het",
  ui     = scBubbHeat_ui,
  server = scBubbHeat_server
)