# TF regulon subnetwork diagram, built from the regulon-to-target-gene table
# id     = "regulon_network"
# title  = "TF Regulon Network"

############################################### Functions ############################################

# subset the regulon table to the selected TF(s)/regulon(s) and build a tbl_graph
# marks nodes belonging to the selected TF set so they can be styled as labels
# rather than plain text in the plot
sc_regulon_build_graph <- function(TF_regulation, tf_sel) {
  edges <- TF_regulation[TF_regulation$group %in% tf_sel, ]
  validate(need(nrow(edges) > 0, "No regulon edges found for the selected TF(s)"))

  tidygraph::as_tbl_graph(edges, directed = FALSE) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(is_tf = name %in% tf_sel)
}

# draws the circular/tree regulon diagram for a given graph
# col_pos/col_neg are hex strings for the positive/negative regulation edges
sc_regulon_plot <- function(g, layout, circular, node_sz, edge_alpha, col_pos, col_neg) {
  ggraph::ggraph(g, layout = layout, circular = circular) +
    ggraph::geom_edge_diagonal(ggplot2::aes(color = regulation), alpha = edge_alpha, width = 1.5) +
    ggraph::geom_node_point(size = node_sz, shape = 21, stroke = 0.5, fill = "grey") +
    ggraph::geom_node_text(ggplot2::aes(label = ifelse(is_tf, "", name)), repel = TRUE, size = 4) +
    ggraph::geom_node_label(ggplot2::aes(label = ifelse(is_tf, name, "")), size = 5,
                            label.padding = ggplot2::unit(0.05, "cm"), label.size = 0.1, repel = TRUE) +
    ggraph::scale_edge_color_manual(values = c("positive" = col_pos, "negative" = col_neg)) +
    ggplot2::scale_x_continuous(expand = c(0.1, 0)) +
    ggplot2::scale_y_continuous(expand = c(0.1, 0), trans = "reverse") +
    ggplot2::theme_void() +
    ggplot2::coord_flip() +
    ggplot2::theme(
      legend.text  = ggplot2::element_text(size = sList["lgd.txt"]),
      legend.title = ggplot2::element_text(size = sList["title"])
    )
}

############################################### UI ###################################################

regulon_network_ui <- function(id, sc1conf, sc1def, ...) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectizeInput(
        inputId  = ns("inpTF"),
        label    = "TF(s) / regulon(s)",
        choices  = NULL,
        selected = NULL,
        multiple = TRUE,
        options  = list(placeholder = "Search TF…")
      ),

      selectInput(
        inputId  = ns("inpLayout"),
        label    = "Layout",
        choices  = c("tree", "kk", "fr", "circle"),
        selected = "tree"
      ),
      checkboxInput(ns("inpCircular"), "Circular layout", value = TRUE),

      sliderInput(ns("inpNodeSz"),    "Node size",    min = 1,   max = 8, value = 3,   step = 0.5),
      sliderInput(ns("inpEdgeAlpha"), "Edge opacity", min = 0.1, max = 1, value = 0.5, step = 0.1),

      fluidRow(
        column(6, textInput(ns("inpColPos"), "Positive edge colour", value = "#4fbbd1")),
        column(6, textInput(ns("inpColNeg"), "Negative edge colour", value = "#b5243e"))
      ),

      hr(),
      fluidRow(
        column(4, selectInput(ns("inpFmt"), "Format", c("png", "pdf", "svg"))),
        column(4, numericInput(ns("inpW"), "Width (in)", 10, min = 1, max = 30)),
        column(4, numericInput(ns("inpH"), "Height (in)", 10, min = 1, max = 30))
      ),
      downloadButton(ns("outDownload"), "Download plot")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Network", plotOutput(ns("outPlot"), height = "700px")),
        tabPanel(
          "Table",
          DT::dataTableOutput(ns("outTable")),
          downloadButton(ns("outTableDL"), "Download table")
        )
      )
    )
  )
}

############################################### Server ##############################################

regulon_network_server <- function(id, sc1conf, sc1def, sc1meta, sc1gene, inpH5, dir_inputs, ...) {
  moduleServer(id, function(input, output, session) {

    # regulon table is independent of prepShinyCellModular — dropped directly
    # into dir_inputs (regulon/), not written by prep. small and static for
    # the session, so read once here and never re-read inside a reactive
    regulon_path <- file.path(dir_inputs, "TF_regulation.rds")
    if (!file.exists(regulon_path)) {
      stop(
        "TF_regulation.rds not found in dir_inputs — place the regulon table (columns: group, item, regulation) there first",
        call. = FALSE
      )
    }
    TF_regulation <- readRDS(regulon_path)

    updateSelectizeInput(
      session,
      inputId  = "inpTF",
      choices  = sort(unique(TF_regulation$group)),
      selected = sort(unique(TF_regulation$group))[1],
      server   = TRUE
    )

    make_plot <- function() {
      req(input$inpTF)
      g <- sc_regulon_build_graph(TF_regulation, input$inpTF)
      sc_regulon_plot(
        g, layout = input$inpLayout, circular = input$inpCircular,
        node_sz = input$inpNodeSz, edge_alpha = input$inpEdgeAlpha,
        col_pos = input$inpColPos, col_neg = input$inpColNeg
      )
    }

    output$outPlot <- renderPlot({
      make_plot()
    })

    output$outDownload <- downloadHandler(
      filename = function() paste0("regulon_network_", Sys.Date(), ".", input$inpFmt),
      content = function(file) {
        # rebuild the plot here — never pull the cached renderPlot object
        p <- make_plot()

        fmt <- tolower(input$inpFmt)
        w   <- input$inpW
        h   <- input$inpH

        if (fmt == "png") {
          png(file, width = w, height = h, units = "in", res = 150)
        } else if (fmt == "pdf") {
          pdf(file, width = w, height = h)
        } else if (fmt == "svg") {
          svg(file, width = w, height = h)
        }
        print(p)
        dev.off()
      }
    )

    tbl_data <- reactive({
      req(input$inpTF)
      TF_regulation[TF_regulation$group %in% input$inpTF, ]
    })

    output$outTable <- DT::renderDataTable({
      DT::datatable(
        tbl_data(),
        rownames = FALSE,
        filter   = "top",
        options  = list(pageLength = 15, scrollX = TRUE)
      )
    })

    output$outTableDL <- downloadHandler(
      filename = function() paste0("regulon_network_table_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(tbl_data(), file, row.names = FALSE)
    )
  })
}

############################################### Registration #########################################

register_tab(
  id          = "regulon_network",
  title       = "TF Regulon Network",
  ui          = regulon_network_ui,
  server      = regulon_network_server,
  author      = "Laura Perlaza-Jimenez",
  description = "Regulon-to-target-gene subnetwork diagram for any selected TF(s), built with tidygraph/ggraph directly from the regulon table",
  version     = "1.0",
  date        = "Aug 2026",
  source      = "TF Regulatory Subnetwork Diagrams chapter",
  contact     = "laura.perlaza-jimenez@monash.edu"
)
