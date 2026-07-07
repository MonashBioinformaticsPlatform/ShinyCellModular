# id     = "genexpr_waffle",
# title  = "Gene Expression Waffle Plot",

############################################### Functions ############################################

# Get gene list, same parsing logic as used elsewhere in the app (kept local here
# so this module stays self-contained and does not depend on other tab files)
scGeneList <- function(inp, inpGene) {

  if (is.null(inp) || is.na(inp) || !nzchar(inp)) {
    return(data.table::data.table(gene = character(0), present = logical(0)))
  }

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

# Bin cells into up to 100 tiles per gene x group, ordered by expression (high to low),
# and average the expression within each tile, adapted from scRNAseqApp (jianhong/scRNAseqApp)
# scWafflePlot(). Groups share the same bin scale (relative to the largest group for that
# gene) so a smaller group naturally fills fewer tiles instead of being stretched to 100.
scWaffleData <- function(ggData) {

  dt <- ggData[, c("geneName", "grpBy", "val"), with = FALSE]
  dt <- dt[!is.na(val)]
  dt <- dt[order(geneName, grpBy, -val)]

  groupSizes <- dt[, .N, by = c("geneName", "grpBy")]
  groupMax   <- groupSizes[, .(N = max(N)), by = "geneName"]
  dt <- merge(dt, groupMax, by = "geneName")

  dt[, rowInGrp := seq_len(.N), by = c("geneName", "grpBy")]
  dt[, idx := cut(seq_len(N[1]), 100, labels = FALSE)[rowInGrp], by = c("geneName", "grpBy")]

  waffleData <- dt[, .(exprs = mean(val, na.rm = TRUE)), by = c("geneName", "grpBy", "idx")]
  waffleData[, x := (idx - 1) %% 10 + 1]
  waffleData[, y := (idx - 1) %/% 10 + 1]
  waffleData[]
}

# Plot gene expression waffle grid, alternative to the bubbleplot / heatmap for a
# handful of genes, adapted from scRNAseqApp (jianhong/scRNAseqApp) scWafflePlot / scDRwafflePlot
scWaffle <- function(inpConf, inpMeta, inp, inpGrp,
                     inpsub1, inpsub2, inpH5, inpGene,
                     inpcols, inpfsz) {

  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]

  geneList <- scGeneList(inp, inpGene)
  geneList <- geneList[present == TRUE]

  shiny::validate(shiny::need(nrow(geneList) >= 1, "Please input at least 1 gene to plot."))
  shiny::validate(shiny::need(nrow(geneList) <= 12,
                              "More than 12 genes to plot. Please reduce the gene list."))

  h5file <- H5File$new(inpH5, mode = "r")
  on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
  h5data <- h5file[["grp"]][["data"]]

  ggData <- data.table::data.table()
  for (iGene in geneList$gene) {
    tmp <- inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) <- c("sampleID", "sub")
    tmp$grpBy <- inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName <- iGene
    tmp$val <- h5data$read(args = list(inpGene[iGene], quote(expr=)))
    ggData <- data.table::rbindlist(list(ggData, tmp))
  }

  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2]
  }

  ggData[val < 0, val := 0]

  shiny::validate(shiny::need(data.table::uniqueN(ggData$grpBy) > 0,
                              "No groups left to plot, please review the filters."))

  waffleData <- scWaffleData(ggData)
  waffleData$geneName <- factor(waffleData$geneName, levels = geneList$gene)

  ggOut <- ggplot(waffleData, aes(x, y, fill = exprs)) +
    geom_tile(color = "white") +
    coord_equal() +
    facet_grid(geneName ~ grpBy) +
    scale_fill_gradientn(colours = cList[[inpcols]]) +
    xlab("") + ylab("") +
    sctheme(base_size = sList[inpfsz], XYval = FALSE) +
    theme(
      panel.grid   = element_blank(),
      strip.text.y = element_text(angle = 0)
    ) +
    guides(fill = guide_colorbar(barwidth = grid::unit(3, "mm")))

  ggOut
}

############################################### UI ###################################################

scWaffle_ui <- function(id, sc1conf, sc1def) {

  ns <- NS(id)

  tabPanel(
    HTML("Waffle Plot"),
    h4("Gene expression waffle plot"),
    "In this tab, users can visualise the expression of a handful of genes across ",
    "groups of cells (e.g. library / cluster) as a grid of unit tiles instead of a ",
    "bubbleplot / heatmap. Each gene x group panel is filled proportionally to how ",
    "many cells that group has, so group size differences are visible at a glance.",
    br(), br(),

    fluidRow(
      column(
        3, style = "border-right: 2px solid black",

        textAreaInput(
          ns("sc1wf1inp"),
          HTML("List of gene names <br />
               ( separated <br />
               by , or ; or newline, max 12):"),
          height = "150px",
          value = paste0(sc1def$genes, collapse = ", ")
        ),

        uiOutput(ns("sc1wf1oupTxt")),
        br(),

        selectInput(
          ns("sc1wf1grp"), "Group by:",
          choices = sc1conf[grp == TRUE]$UI,
          selected = sc1conf[grp == TRUE]$UI[1]
        ),

        radioButtons(
          ns("sc1wf1cols"), "Colour scheme:",
          choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
          selected = "White-Red"
        ),

        br(),

        actionButton(ns("sc1wf1togL"), "Filter Cells"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1wf1togL")),
          selectInput(ns("sc1wf1sub1"), "Cell information to subset:",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp1),
          uiOutput(ns("sc1wf1sub1.ui")),
          actionButton(ns("sc1wf1sub1all"), "Select all groups", class = "btn btn-primary"),
          actionButton(ns("sc1wf1sub1non"), "Deselect all groups", class = "btn btn-primary")
        ),

        br(), br(),

        actionButton(ns("sc1wf1tog"), "Customize Plot"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1wf1tog")),
          radioButtons(ns("sc1wf1psz"), "Plot size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Medium", inline = TRUE),
          radioButtons(ns("sc1wf1fsz"), "Font size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Medium", inline = TRUE)
        )
      ),

      column(
        9,
        uiOutput(ns("sc1wf1oup.ui")),
        downloadButton(ns("sc1wf1oup.pdf"), "Download PDF"),
        downloadButton(ns("sc1wf1oup.png"), "Download PNG"),
        br(),

        div(style = "display:inline-block",
            numericInput(ns("sc1wf1oup.h"), "PDF / PNG height:", width = "138px",
                         min = 4, max = 20, value = 8, step = 0.5)),

        div(style = "display:inline-block",
            numericInput(ns("sc1wf1oup.w"), "PDF / PNG width:", width = "138px",
                         min = 4, max = 20, value = 10, step = 0.5))
      )
    )
  )
}

############################################### Server ##############################################

scWaffle_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    if (!exists("pList2", inherits = TRUE)) {
      pList2 <<- c(Small = "450px", Medium = "650px", Large = "850px")
    }

    # Output text, same feedback style as the bubbleplot / heatmap tab
    output$sc1wf1oupTxt <- renderUI({
      geneList <- scGeneList(input$sc1wf1inp, sc1gene)

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

    # Filter UI (independent)
    output$sc1wf1sub1.ui <- renderUI({
      req(input$sc1wf1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1wf1sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(ns("sc1wf1sub2"), "Select which cells to show",
                         inline = TRUE, choices = sub, selected = sub)
    })

    observeEvent(input$sc1wf1sub1non, {
      req(input$sc1wf1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1wf1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1wf1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = NULL, inline = TRUE)
    })

    observeEvent(input$sc1wf1sub1all, {
      req(input$sc1wf1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1wf1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1wf1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = sub, inline = TRUE)
    })

    plotFun <- reactive({
      req(input$sc1wf1inp, input$sc1wf1grp)
      scWaffle(
        sc1conf, sc1meta,
        input$sc1wf1inp, input$sc1wf1grp,
        input$sc1wf1sub1, input$sc1wf1sub2,
        file.path(dir_inputs, "RNA", "sc1gexpr.h5"),
        sc1gene,
        input$sc1wf1cols, input$sc1wf1fsz
      )
    })

    output$sc1wf1oup <- renderPlot({
      plotFun()
    })

    output$sc1wf1oup.ui <- renderUI({
      req(input$sc1wf1psz)
      plotOutput(ns("sc1wf1oup"), height = pList2[input$sc1wf1psz])
    })

    output$sc1wf1oup.pdf <- downloadHandler(
      filename = function() {
        paste0("sc1waffle_", input$sc1wf1grp, ".pdf")
      },
      content = function(file) {
        ggsave(
          file, device = "pdf",
          height = input$sc1wf1oup.h, width = input$sc1wf1oup.w,
          useDingbats = FALSE,
          plot = plotFun()
        )
      }
    )

    output$sc1wf1oup.png <- downloadHandler(
      filename = function() {
        paste0("sc1waffle_", input$sc1wf1grp, ".png")
      },
      content = function(file) {
        ggsave(
          file, device = "png",
          height = input$sc1wf1oup.h, width = input$sc1wf1oup.w,
          plot = plotFun()
        )
      }
    )

  })
}

############################################### Registration #########################################

register_tab(
  id          = "genexpr_waffle",
  title       = "Gene Expression Waffle Plot",
  ui          = scWaffle_ui,
  server      = scWaffle_server,
  author      = "Laura Perlaza-Jimenez",
  description = "Gene expression waffle plot, alternative to the bubbleplot / heatmap for a small gene list, adapted from scRNAseqApp (jianhong/scRNAseqApp)",
  version     = "1.0",
  date        = "Jul 2026",
  source      = "scRNAseqApp (jianhong/scRNAseqApp), MGBP custom",
  contact     = "laura.perlaza-jimenez@monash.edu"
)
