# Trajectory / pseudotime plot tab
# id     = "trajectory"
# title  = "Trajectory"
# requires: trajectory_methods.rds and trajectory_<method>.rds in dir_inputs/
#           (written by prepShinyCellModular() when do_trajectory = TRUE)
# pseudotime values are read directly from sc1meta as pseudotime_<method>

############################################### Functions ############################################

sc_traj_available_methods <- function(dir_inputs) {
  manifest_file <- file.path(dir_inputs, "trajectory_methods.rds")
  if (!file.exists(manifest_file)) return(character(0))
  readRDS(manifest_file)
}

sc_traj_load_graph <- function(dir_inputs, method) {
  graph_file <- file.path(dir_inputs, paste0("trajectory_", method, ".rds"))
  if (!file.exists(graph_file)) return(NULL)
  readRDS(graph_file)
}

# 2D pseudotime plot with trajectory graph overlay, optionally split by a group column
sc_traj_plot2d <- function(sc1conf, sc1meta, dir_inputs, method,
                           inpdrX, inpdrY, inpsplit, inpsiz, inpcol, inpfsz, inptxt) {

  entry <- sc_traj_load_graph(dir_inputs, method)
  shiny::validate(shiny::need(!is.null(entry), paste0("No trajectory graph found for method: ", method)))

  pt_col <- paste0("pseudotime_", method)
  shiny::validate(shiny::need(pt_col %in% names(sc1meta),
                               paste0("Column ", pt_col, " not found in metadata.")))

  x_id <- sc1conf[UI == inpdrX]$ID
  y_id <- sc1conf[UI == inpdrY]$ID

  cols <- c(x_id, y_id, pt_col)
  split_id <- NULL
  if (!is.null(inpsplit) && nzchar(inpsplit) && inpsplit != "None") {
    split_id <- sc1conf[UI == inpsplit]$ID
    cols <- c(cols, split_id)
  }

  ggData <- sc1meta[, cols, with = FALSE]
  data.table::setnames(ggData, c("X", "Y", "pseudotime")[seq_len(3)])
  if (!is.null(split_id)) data.table::setnames(ggData, ncol(ggData), "split")

  p <- ggplot2::ggplot(ggData, ggplot2::aes(X, Y, color = pseudotime)) +
    ggplot2::geom_point(size = inpsiz, shape = 16) +
    ggplot2::geom_segment(
      data = entry$graph,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      inherit.aes = FALSE, color = "black", linewidth = 0.4
    ) +
    ggplot2::scale_color_gradientn("Pseudotime", colours = cList[[inpcol]]) +
    ggplot2::xlab(inpdrX) + ggplot2::ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt)

  if (!is.null(split_id)) {
    p <- p + ggplot2::facet_wrap(~split)
  }

  p
}

# 3D pseudotime plot with trajectory graph overlay, reuses the scGet3DReduction
# helper already defined in cellinfo3D_cellinfo3D.R (loaded earlier at app build time)
sc_traj_plot3d <- function(sc1conf, sc1meta, dir_inputs, method, dr3d_base, inpsiz, inpcol) {

  entry <- sc_traj_load_graph(dir_inputs, method)
  shiny::validate(shiny::need(!is.null(entry), paste0("No trajectory graph found for method: ", method)))

  pt_col <- paste0("pseudotime_", method)
  shiny::validate(shiny::need(pt_col %in% names(sc1meta),
                               paste0("Column ", pt_col, " not found in metadata.")))

  red3 <- scGet3DReduction(sc1conf, sc1meta, base = dr3d_base)

  ggData <- sc1meta[, c(red3$ID$X, red3$ID$Y, red3$ID$Z, pt_col), with = FALSE]
  data.table::setnames(ggData, c("X", "Y", "Z", "pseudotime"))

  p <- plotly::plot_ly()

  p <- p %>% plotly::add_trace(
    data = ggData, x = ~X, y = ~Y, z = ~Z,
    type = "scatter3d", mode = "markers",
    marker = list(size = inpsiz, color = ~pseudotime,
                  colorscale = cList[[inpcol]], showscale = TRUE),
    hovertemplate = paste0("Pseudotime: %{marker.color:.2f}<extra></extra>")
  )

  # graph overlay: one line trace per segment (kept simple and robust;
  # segments are typically few enough that per-row traces are not a
  # performance concern, unlike per-cell markers)
  g <- entry$graph
  for (i in seq_len(nrow(g))) {
    p <- p %>% plotly::add_trace(
      x = c(g$x_start[i], g$x_end[i]),
      y = c(g$y_start[i], g$y_end[i]),
      z = c(g$z_start[i], g$z_end[i]),
      type = "scatter3d", mode = "lines",
      line = list(color = "black", width = 4),
      showlegend = FALSE, hoverinfo = "skip"
    )
  }

  p <- p %>%
    plotly::layout(
      scene = list(
        xaxis = list(title = red3$UI$X),
        yaxis = list(title = red3$UI$Y),
        zaxis = list(title = red3$UI$Z)
      )
    ) %>%
    plotly::config(toImageButtonOptions = list(
      format = "svg", filename = "trajectory_3d", width = 800, height = 600
    ))

  list(plot = p, reduction = red3$base)
}

############################################### UI ###################################################

trajectory_ui <- function(id, sc1conf, sc1def) {

  ns <- NS(id)

  tabPanel(
    HTML("Trajectory"),
    h4("Pseudotime trajectory"),
    "Colours cells by pseudotime and overlays the fitted trajectory graph. ",
    "Requires trajectory data written by prepShinyCellModular(do_trajectory = TRUE).",
    br(), br(),

    fluidRow(
      column(
        3, style = "border-right: 2px solid black",

        radioButtons(ns("sc1traj_mode"), "Dimensions:",
                     choices = c("2D", "3D"), selected = "2D", inline = TRUE),

        uiOutput(ns("sc1traj_method_ui")),

        conditionalPanel(
          condition = sprintf("input['%s'] == '2D'", ns("sc1traj_mode")),
          selectInput(ns("sc1traj_dimX"), "X-axis:", choices = sc1conf[dimred == TRUE]$UI),
          selectInput(ns("sc1traj_dimY"), "Y-axis:", choices = sc1conf[dimred == TRUE]$UI),
          selectInput(ns("sc1traj_split"), "Split by (optional):",
                      choices = c("None", sc1conf[grp == TRUE]$UI), selected = "None")
        ),

        conditionalPanel(
          condition = sprintf("input['%s'] == '3D'", ns("sc1traj_mode")),
          uiOutput(ns("sc1traj_3d_reduction_ui"))
        ),

        br(),
        sliderInput(ns("sc1traj_siz"), "Point size:", min = 0.1, max = 5, value = 1.25, step = 0.1),
        selectInput(ns("sc1traj_col"), "Colour scheme:", choices = names(cList), selected = names(cList)[1]),
        radioButtons(ns("sc1traj_fsz"), "Font size:",
                     choices = c("Small", "Medium", "Large"), selected = "Medium", inline = TRUE),
        checkboxInput(ns("sc1traj_txt"), "Show axis text", value = TRUE)
      ),

      column(
        9,
        conditionalPanel(
          condition = sprintf("input['%s'] == '2D'", ns("sc1traj_mode")),
          plotOutput(ns("sc1traj_plot2d"), height = "600px")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == '3D'", ns("sc1traj_mode")),
          plotly::plotlyOutput(ns("sc1traj_plot3d"), height = "600px")
        ),
        br(),
        div(
          style = "display:flex; gap:8px; align-items:center;",
          downloadButton(ns("sc1traj_pdf"), "PDF"),
          downloadButton(ns("sc1traj_png"), "PNG")
        )
      )
    )
  )
}

############################################### Server ##############################################

trajectory_server <- function(id, sc1conf, sc1def, sc1meta, dir_inputs) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    available_methods <- sc_traj_available_methods(dir_inputs)

    output$sc1traj_method_ui <- renderUI({
      shiny::validate(shiny::need(length(available_methods) > 0,
                                   "No trajectory data found. Run prepShinyCellModular(do_trajectory = TRUE) with a trajectory_* entry in @misc first."))
      selectInput(ns("sc1traj_method"), "Trajectory method:",
                  choices = available_methods, selected = available_methods[1])
    })

    output$sc1traj_3d_reduction_ui <- renderUI({
      dt <- scList3DReductions(sc1conf, sc1meta)
      shiny::validate(shiny::need(nrow(dt) > 0, "No 3D reduction available. Run prepShinyCellModular(do_umap3d = TRUE) first."))
      selectInput(ns("sc1traj_3d_reduction"), "3D reduction:", choices = dt$base, selected = dt$base[1])
    })

    plot2d <- reactive({
      req(input$sc1traj_method, input$sc1traj_dimX, input$sc1traj_dimY)
      sc_traj_plot2d(
        sc1conf, sc1meta, dir_inputs, input$sc1traj_method,
        input$sc1traj_dimX, input$sc1traj_dimY, input$sc1traj_split,
        input$sc1traj_siz, input$sc1traj_col, input$sc1traj_fsz, input$sc1traj_txt
      )
    })

    plot3d <- reactive({
      req(input$sc1traj_method, input$sc1traj_3d_reduction)
      sc_traj_plot3d(
        sc1conf, sc1meta, dir_inputs, input$sc1traj_method,
        input$sc1traj_3d_reduction, input$sc1traj_siz, input$sc1traj_col
      )
    })

    output$sc1traj_plot2d <- renderPlot({ plot2d() })
    output$sc1traj_plot3d <- plotly::renderPlotly({ plot3d()$plot })

    output$sc1traj_pdf <- downloadHandler(
      filename = function() paste0("trajectory_", input$sc1traj_method, ".pdf"),
      content = function(file) {
        req(input$sc1traj_mode == "2D")
        ggplot2::ggsave(file, plot = plot2d(), device = "pdf", width = 7, height = 6, useDingbats = FALSE)
      }
    )

    output$sc1traj_png <- downloadHandler(
      filename = function() paste0("trajectory_", input$sc1traj_method, ".png"),
      content = function(file) {
        req(input$sc1traj_mode == "2D")
        ggplot2::ggsave(file, plot = plot2d(), device = "png", width = 7, height = 6, dpi = 300)
      }
    )
  })
}

############################################### Registration #########################################

register_tab(
  id          = "trajectory_multi",
  title       = "Trajectory (Multi)",
  ui          = trajectory_ui,
  server      = trajectory_server,
  author      = "Laura Perlaza-Jimenez",
  description = "Pseudotime trajectory plot with Monocle3 or Slingshot graph overlay, new tab created by MGBP",
  version     = "1.0",
  date        = "Jul 2026",
  source      = "MGBP custom",
  contact     = "laura.perlaza-jimenez@monash.edu"
)
