# id     = "split_violin",
# title  = "Split Violin Plot",

############################################### Functions ############################################

# Split violin geom, lets 2 groups be compared side by side within one violin shape
# instead of 2 separate violins, adapted from scRNAseqApp (jianhong/scRNAseqApp)
# original source: jan-glx on StackOverflow
# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data$xminv <- data$x - data$violinwidth * (data$x - data$xmin)
    data$xmaxv <- data$x + data$violinwidth * (data$xmax - data$x)
    grp <- data[1, 'group']
    if (grp %% 2 == 1) {
      data$x <- data$xminv
      data.order <- data$y
    } else {
      data$x <- data$xmaxv
      data.order <- -data$y
    }
    newdata <- data[order(data.order), , drop = FALSE]
    newdata <- rbind(
      newdata[1, ],
      newdata,
      newdata[nrow(newdata), ],
      newdata[1, ]
    )
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), 'x'] <- round(newdata[1, 'x'])
    grob <- if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::QuantileSegments(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile.grob <- GeomPath$draw_panel(both, ...)
      grid::grobTree(GeomPolygon$draw_panel(newdata, ...), name = quantile.grob)
    } else {
      GeomPolygon$draw_panel(newdata, ...)
    }
    grob$name <- grid::grobName(grob, prefix = 'geom_split_violin')
    grob
  }
)

# layer constructor for the split violin geom, same usage as geom_violin()
geom_split_violin <- function(mapping = NULL, data = NULL, stat = 'ydensity',
                              position = 'identity', ..., draw_quantiles = NULL,
                              trim = TRUE, scale = 'area', na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...)
  )
}

# Plot gene expression / cell info as a split violin, compares 2 groups per X category
scSplitViolin <- function(inpConf, inpMeta, inp1, inp2, inp3,
                          inpLvlA, inpLvlB,
                          inpsub1, inpsub2, inpH5, inpGene,
                          inppts, inpsiz, inpfsz, inscale_min, inscale_max,
                          x_order = NULL) {

  shiny::validate(shiny::need(!is.null(inpLvlA) && !is.null(inpLvlB) && inpLvlA != inpLvlB,
                              "Please choose 2 different groups to split by."))

  expr_min <- inscale_min
  expr_max <- inscale_max

  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]

  ggData <- inpMeta[, c(
    inpConf[UI == inp1]$ID,
    inpConf[UI == inp2]$ID,
    inpConf[UI == inpsub1]$ID
  ), with = FALSE]
  colnames(ggData) <- c("X", "splitGrp", "sub")

  # keep only the 2 chosen split levels, ordered A (left half) then B (right half)
  ggData <- ggData[as.character(splitGrp) %in% c(inpLvlA, inpLvlB)]
  ggData$splitGrp <- factor(as.character(ggData$splitGrp), levels = c(inpLvlA, inpLvlB))

  if (inp3 %in% inpConf$UI) {
    ggData$val <- inpMeta[[inpConf[UI == inp3]$ID]]
  } else {
    h5file <- H5File$new(inpH5, mode = "r")
    on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
    h5data <- h5file[["grp"]][["data"]]
    ggData$val <- h5data$read(args = list(inpGene[inp3], quote(expr=)))
    ggData[val < 0]$val <- 0
    set.seed(42)
    tmpNoise <- rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000
    ggData$val <- ggData$val + tmpNoise
  }

  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2]
  }

  shiny::validate(shiny::need(nrow(ggData) > 0, "No cells left to plot, please review the filters."))

  ggLvl <- levels(inpMeta[[inpConf[UI == inp1]$ID]])
  if (is.null(ggLvl)) ggLvl <- sort(unique(as.character(ggData$X)))
  ggLvl <- ggLvl[ggLvl %in% unique(as.character(ggData$X))]

  if (!is.null(x_order) && length(x_order) > 0) {
    ggLvl <- x_order[x_order %in% ggLvl]
  }

  ggData$X <- factor(as.character(ggData$X), levels = ggLvl)

  ggOut <- ggplot(ggData, aes(X, val, fill = splitGrp)) +
    geom_split_violin(scale = "width", trim = TRUE)

  if (isTRUE(inppts)) {
    ggOut <- ggOut + geom_jitter(size = inpsiz, shape = 16, width = 0.1, alpha = 0.5)
  }

  ggOut <- ggOut +
    xlab(inp1) + ylab(inp3) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_fill_manual(paste0(inp2, ":"), values = c("#5773CC", "#FFB900")) +
    theme(legend.position = "right")

  if (!is.null(expr_min) && !is.null(expr_max) &&
      !is.na(expr_min) && !is.na(expr_max)) {
    ggOut <- ggOut + scale_y_continuous(limits = c(expr_min, expr_max))
  }

  ggOut
}

############################################### UI ###################################################

scSplitViolin_ui <- function(id, sc1conf, sc1def) {

  ns <- NS(id)

  tabPanel(
    HTML("Split Violin Plot"),
    h4("Compare 2 groups side by side within one violin"),
    "In this tab, users can visualise the gene expression or continuous cell information ",
    "(Y-axis) across groups of cells (X-axis), split by a second 2-level grouping (e.g. ",
    "treatment vs control), shown as left/right halves of the same violin.",
    br(), br(),

    fluidRow(
      column(
        3, style = "border-right: 2px solid black",

        selectInput(
          ns("sc1sv1inp1"), "Cell information (X-axis):",
          choices = sc1conf[grp == TRUE]$UI,
          selected = sc1def$grp1
        ),

        selectInput(
          ns("sc1sv1inp2"), "Split by (2 groups shown as left/right half):",
          choices = sc1conf[grp == TRUE]$UI,
          selected = sc1conf[grp == TRUE]$UI[min(2, length(sc1conf[grp == TRUE]$UI))]
        ),

        uiOutput(ns("sc1sv1splitLvl.ui")),

        selectInput(ns("sc1sv1inp3"), "Cell Info / Gene name (Y-axis):", choices = NULL),

        checkboxInput(ns("sc1sv1pts"), "Show data points", value = FALSE),

        actionButton(ns("sc1sv1togOrderX"), "Order X axis"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1sv1togOrderX")),
          h5("Drag to reorder X axis groups"),
          uiOutput(ns("sc1sv1xorder.ui")),
          actionButton(ns("sc1sv1xorder_reset"), "Reset to default", class = "btn btn-primary")
        ),

        br(), br(),

        actionButton(ns("sc1sv1togL"), "Filter Cells"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1sv1togL")),
          selectInput(ns("sc1sv1sub1"), "Cell information to subset:",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp1),
          uiOutput(ns("sc1sv1sub1.ui")),
          actionButton(ns("sc1sv1sub1all"), "Select all groups", class = "btn btn-primary"),
          actionButton(ns("sc1sv1sub1non"), "Deselect all groups", class = "btn btn-primary")
        ),

        br(), br(),

        actionButton(ns("sc1sv1tog"), "Customize Plot"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1sv1tog")),
          sliderInput(ns("sc1sv1siz"), "Data point size:",
                      min = 0, max = 4, value = 1.25, step = 0.25),
          radioButtons(ns("sc1sv1psz"), "Plot size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Medium", inline = TRUE),
          radioButtons(ns("sc1sv1fsz"), "Font size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Medium", inline = TRUE)
        ),

        actionButton(ns("sc1sv1fixscale"), "Fix Y scale", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("sc1sv1fixscale")),
          fluidRow(
            column(6, numericInput(ns("sc1sv1ymin"), "Y min", value = NULL, step = 0.1)),
            column(6, numericInput(ns("sc1sv1ymax"), "Y max", value = NULL, step = 0.1))
          )
        )
      ),

      column(
        9,
        uiOutput(ns("sc1sv1oup.ui")),
        downloadButton(ns("sc1sv1oup.pdf"), "Download PDF"),
        downloadButton(ns("sc1sv1oup.png"), "Download PNG"),
        br(),

        div(style = "display:inline-block",
            numericInput(ns("sc1sv1oup.h"), "PDF / PNG height:", width = "138px",
                         min = 4, max = 20, value = 8, step = 0.5)),

        div(style = "display:inline-block",
            numericInput(ns("sc1sv1oup.w"), "PDF / PNG width:", width = "138px",
                         min = 4, max = 20, value = 10, step = 0.5))
      )
    )
  )
}

############################################### Server ##############################################

scSplitViolin_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns
    observe_helpers()

    optCrt <- "{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"

    updateSelectizeInput(
      session, "sc1sv1inp3", server = TRUE,
      choices = c(sc1conf[is.na(fID)]$UI, names(sc1gene)),
      selected = sc1conf[is.na(fID)]$UI[1],
      options = list(
        maxOptions = length(sc1conf[is.na(fID)]$UI) + 3,
        create = TRUE, persist = TRUE, render = I(optCrt)
      )
    )

    if (!exists("pList2", inherits = TRUE)) {
      pList2 <<- c(Small = "450px", Medium = "650px", Large = "850px")
    }

    # Split levels UI, depends on which cell info is picked to split by
    output$sc1sv1splitLvl.ui <- renderUI({
      req(input$sc1sv1inp2)
      lvl <- strsplit(sc1conf[UI == input$sc1sv1inp2]$fID, "\\|")[[1]]
      tagList(
        selectInput(ns("sc1sv1lvlA"), "Left half:", choices = lvl, selected = lvl[1]),
        selectInput(ns("sc1sv1lvlB"), "Right half:", choices = lvl,
                    selected = lvl[min(2, length(lvl))])
      )
    })

    # Filter UI (independent)
    output$sc1sv1sub1.ui <- renderUI({
      req(input$sc1sv1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1sv1sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(ns("sc1sv1sub2"), "Select which cells to show",
                         inline = TRUE, choices = sub, selected = sub)
    })

    observeEvent(input$sc1sv1sub1non, {
      req(input$sc1sv1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1sv1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1sv1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = NULL, inline = TRUE)
    })

    observeEvent(input$sc1sv1sub1all, {
      req(input$sc1sv1sub1)
      sub <- strsplit(sc1conf[UI == input$sc1sv1sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "sc1sv1sub2",
                               label = "Select which cells to show",
                               choices = sub, selected = sub, inline = TRUE)
    })

    # X axis ordering levels come from the X axis variable itself (independent from filter)
    x_levels_default <- reactive({
      req(input$sc1sv1inp1)
      x_id <- sc1conf[UI == input$sc1sv1inp1]$ID
      req(length(x_id) == 1, !is.na(x_id), x_id != "")
      levs <- levels(sc1meta[[x_id]])
      if (is.null(levs)) levs <- sort(unique(as.character(sc1meta[[x_id]])))
      levs
    })

    x_order <- reactiveVal(NULL)

    observeEvent(x_levels_default(), {
      x_order(x_levels_default())
    }, ignoreInit = FALSE)

    output$sc1sv1xorder.ui <- renderUI({
      req(x_levels_default())
      ord <- x_order()
      if (is.null(ord) || length(ord) == 0) ord <- x_levels_default()
      ord <- ord[ord %in% x_levels_default()]
      if (length(ord) == 0) ord <- x_levels_default()

      sortable::rank_list(
        text = "X axis group order",
        labels = ord,
        input_id = ns("sc1sv1xorder_rank")
      )
    })

    observeEvent(input$sc1sv1xorder_rank, {
      req(input$sc1sv1xorder_rank)
      x_order(input$sc1sv1xorder_rank)
    })

    observeEvent(input$sc1sv1xorder_reset, {
      x_order(x_levels_default())
    })

    x_order_final <- reactive({
      levs <- x_levels_default()
      ord <- x_order()
      if (is.null(ord) || length(ord) == 0) return(levs)
      ord <- ord[ord %in% levs]
      if (length(ord) == 0) levs else ord
    })

    plotFun <- reactive({
      req(input$sc1sv1inp1, input$sc1sv1inp2, input$sc1sv1inp3,
          input$sc1sv1lvlA, input$sc1sv1lvlB)
      scSplitViolin(
        sc1conf, sc1meta,
        input$sc1sv1inp1, input$sc1sv1inp2, input$sc1sv1inp3,
        input$sc1sv1lvlA, input$sc1sv1lvlB,
        input$sc1sv1sub1, input$sc1sv1sub2,
        file.path(dir_inputs, "RNA", "sc1gexpr.h5"),
        sc1gene,
        input$sc1sv1pts, input$sc1sv1siz, input$sc1sv1fsz,
        input$sc1sv1ymin, input$sc1sv1ymax,
        x_order = x_order_final()
      )
    })

    output$sc1sv1oup <- renderPlot({
      plotFun()
    })

    output$sc1sv1oup.ui <- renderUI({
      req(input$sc1sv1psz)
      plotOutput(ns("sc1sv1oup"), height = pList2[input$sc1sv1psz])
    })

    output$sc1sv1oup.pdf <- downloadHandler(
      filename = function() {
        paste0("sc1splitviolin_", input$sc1sv1inp1, "_", input$sc1sv1inp3, ".pdf")
      },
      content = function(file) {
        ggsave(
          file, device = "pdf",
          height = input$sc1sv1oup.h, width = input$sc1sv1oup.w,
          useDingbats = FALSE,
          plot = plotFun()
        )
      }
    )

    output$sc1sv1oup.png <- downloadHandler(
      filename = function() {
        paste0("sc1splitviolin_", input$sc1sv1inp1, "_", input$sc1sv1inp3, ".png")
      },
      content = function(file) {
        ggsave(
          file, device = "png",
          height = input$sc1sv1oup.h, width = input$sc1sv1oup.w,
          plot = plotFun()
        )
      }
    )

  })
}

############################################### Registration #########################################

register_tab(
  id          = "split_violin",
  title       = "Split Violin Plot",
  ui          = scSplitViolin_ui,
  server      = scSplitViolin_server,
  author      = "Laura Perlaza-Jimenez",
  description = "Split violin plot comparing 2 groups side by side within one violin, adapted from scRNAseqApp (jianhong/scRNAseqApp), original geom by jan-glx (StackOverflow)",
  version     = "1.0",
  date        = "Jul 2026",
  source      = "scRNAseqApp (jianhong/scRNAseqApp), MGBP custom",
  contact     = "laura.perlaza-jimenez@monash.edu"
)
