# CellInfo3D vs CellInfo3D
# This tab is base in the orginal "cellinfo_cellinfo" in ShinyCell
# id     = "cellinfo3D_cellinfo3D",
# title  = "CellInfo3D vs CellInfo3D",

########################################## Functions ########################################################

scList3DReductions <- function(inpConf, inpMeta) {
  
  conf <- data.table::as.data.table(inpConf)
  meta_cols <- names(inpMeta)
  
  if (!("dimred" %in% names(conf))) stop("inpConf must contain a 'dimred' column.")
  if (!all(c("UI","ID") %in% names(conf))) stop("inpConf must contain columns UI and ID.")
  
  conf <- conf[dimred == TRUE]
  conf <- conf[!is.na(UI) & !is.na(ID)]
  conf <- conf[UI != "" & ID != ""]
  conf <- conf[ID %in% meta_cols]
  
  if (nrow(conf) == 0) return(data.table::data.table())
  
  conf[, dim := data.table::fifelse(grepl("([._-])?1$", UI), 1L,
                                    data.table::fifelse(grepl("([._-])?2$", UI), 2L,
                                                        data.table::fifelse(grepl("([._-])?3$", UI), 3L, NA_integer_)))]
  
  conf <- conf[!is.na(dim)]
  if (nrow(conf) == 0) return(data.table::data.table())
  
  conf[, base := gsub("([._-])?[123]$", "", UI)]
  
  bases <- conf[, .(has1 = any(dim == 1L), has2 = any(dim == 2L), has3 = any(dim == 3L)),
                by = base][has1 & has2 & has3, .(base)]
  
  if (nrow(bases) == 0) return(data.table::data.table())
  
  out <- conf[base %in% bases$base, .(UI = UI[1], ID = ID[1]), by = .(base, dim)]
  
  out_ui <- data.table::dcast(out, base ~ dim, value.var = "UI")
  out_id <- data.table::dcast(out, base ~ dim, value.var = "ID")
  data.table::setnames(out_ui, c("base", "UI_1", "UI_2", "UI_3"))
  data.table::setnames(out_id, c("base", "ID_1", "ID_2", "ID_3"))
  
  out_ui[out_id, on = "base"][]
}

scGet3DReduction <- function(inpConf, inpMeta, base = NULL) {
  
  dt <- scList3DReductions(inpConf, inpMeta)
  if (nrow(dt) == 0) stop("No 3D reduction found.")
  
  base_sel <- base
  if (is.null(base_sel) || !base_sel %in% dt$base) base_sel <- dt$base[1]
  row <- dt[base == base_sel][1]
  
  list(
    base = row$base,
    UI = list(X = row$UI_1, Y = row$UI_2, Z = row$UI_3),
    ID = list(X = row$ID_1, Y = row$ID_2, Z = row$ID_3)
  )
}

scDRcell3D <- function(inpConf, inpMeta, dr3d_base,
                       inp1, inpsub1, inpsub2,
                       inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab) {
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  red3 <- scGet3DReduction(inpConf, inpMeta, base = dr3d_base)
  
  ggData <- inpMeta[, c(
    red3$ID$X, red3$ID$Y, red3$ID$Z,
    inpConf[UI == inp1]$ID,
    inpConf[UI == inpsub1]$ID
  ), with = FALSE]
  data.table::setnames(ggData, c("X", "Y", "Z", "val", "sub"))
  
  bgCells <- FALSE
  if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData  <- ggData[sub %in% inpsub2]
  }
  
  if (inpord == "Max-1st")       ggData <- ggData[order(val)]
  else if (inpord == "Min-1st")  ggData <- ggData[order(-val)]
  else if (inpord == "Random")   ggData <- ggData[sample(nrow(ggData))]
  
  isDiscrete <- !is.na(inpConf[UI == inp1]$fCL)
  
  if (isDiscrete) {
    ggCol <- strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
    names(ggCol) <- levels(ggData$val)
    ggLvl <- levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)]
    ggData$val <- factor(ggData$val, levels = ggLvl)
    ggCol <- ggCol[ggLvl]
  }
  
  ggData$text <- as.character(ggData$val)
  if (bgCells) ggData2$text <- as.character(ggData2$val)
  
  hovertemplate <- paste(
    "Cluster: %{text}<br>",
    paste0(red3$UI$X, ": %{x:.2f}<br>"),
    paste0(red3$UI$Y, ": %{y:.2f}<br>"),
    paste0(red3$UI$Z, ": %{z:.2f}"),
    "<extra></extra>"
  )
  
  p <- plotly::plot_ly()
  
  if (bgCells) {
    p <- p %>% plotly::add_trace(
      data = ggData2, x = ~X, y = ~Y, z = ~Z,
      type = "scatter3d", mode = "markers",
      marker = list(size = inpsiz, color = "snow2"),
      text = ~text, hovertemplate = hovertemplate, showlegend = FALSE
    )
  }
  
  if (isDiscrete) {
    p <- p %>% plotly::add_trace(
      data = ggData, x = ~X, y = ~Y, z = ~Z,
      type = "scatter3d", mode = "markers",
      color = ~val, colors = unname(ggCol),
      text = ~text, marker = list(size = inpsiz), hovertemplate = hovertemplate
    )
  } else {
    p <- p %>% plotly::add_trace(
      data = ggData, x = ~X, y = ~Y, z = ~Z,
      type = "scatter3d", mode = "markers", text = ~text,
      marker = list(size = inpsiz, color = ggData$val,
                    colorscale = cList[[inpcol]], showscale = TRUE),
      hovertemplate = hovertemplate, showlegend = FALSE
    )
  }
  
  p <- p %>%
    plotly::layout(
      legend = list(itemsizing = "constant", font = list(size = 12)),
      scene = list(
        xaxis = list(title = red3$UI$X),
        yaxis = list(title = red3$UI$Y),
        zaxis = list(title = red3$UI$Z)
      )
    ) %>%
    plotly::config(toImageButtonOptions = list(
      format = "svg", filename = "my_plot", width = 800, height = 600
    ))
  
  list(plot = p, reduction = red3$base)
}

########################################## UI ########################################################

scDRcell3D_ui <- function(id, sc1conf, sc1def) {
  ns <- NS(id)
  
  tabPanel(
    HTML("CellInfo3D vs CellInfo3D"),
    h4("Cell information vs cell information on dimension reduction 3D"),
    "In this tab, users can visualise two cell informations side-by-side ",
    "on low-dimensional representions.",
    br(), br(),
    
    fluidRow(
      column(
        3, h4("Dimension Reduction 3D"),
        selectInput(ns("sc1a2dr3d_base"), "reduction with more than 2 dimensions:",
                    choices = character(0))
      ),
      
      column(
        3,
        # ── Camera sync: Capture & Apply workflow ──
        h5("Sync 3D camera"),
        radioButtons(
          ns("sc1a2master"), "Master plot:",
          choices  = c("Cell Info 1" = "plot1", "Cell Info 2" = "plot2"),
          selected = "plot1", inline = TRUE
        ),
        actionButton(ns("sc1a2captureBtn"), "1. Capture view",
                     class = "btn btn-info", icon = icon("camera")),
        actionButton(ns("sc1a2applyBtn"), "2. Apply to other plot",
                     class = "btn btn-success", icon = icon("sync")),
        br(),
        verbatimTextOutput(ns("sc1a2camTxt")),
        
        br(),
        actionButton(ns("sc1a2togL"), "Filter cells"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2togL")),
          selectInput(ns("sc1a2sub1"), "Cell information to subset:",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp1),
          uiOutput(ns("sc1a2sub1.ui")),
          actionButton(ns("sc1a2sub1all"), "Select all groups", class = "btn btn-primary"),
          actionButton(ns("sc1a2sub1non"), "Deselect all groups", class = "btn btn-primary")
        )
      ),
      
      column(
        6,
        actionButton(ns("sc1a2tog0"), "Customize Aesthetics for Both Plots"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2tog0")),
          fluidRow(
            column(6,
              sliderInput(ns("sc1a2siz"), "Point size:", min = 0, max = 4, value = 1.25, step = 0.25),
              radioButtons(ns("sc1a2psz"), "Plot size:", choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE),
              radioButtons(ns("sc1a2fsz"), "Font size:", choices = c("Small", "Medium", "Large"),
                           selected = "Medium", inline = TRUE)
            ),
            column(6,
              radioButtons(ns("sc1a2asp"), "Aspect ratio:", choices = c("Square", "Fixed", "Free"),
                           selected = "Square", inline = TRUE),
              checkboxInput(ns("sc1a2txt"), "Show axis text", value = FALSE)
            )
          )
        )
      )
    ),
    
    # ── JavaScript: "Capture" reads the LIVE camera from the plotly WebGL scene ──
    tags$script(HTML(sprintf("
      $(document).ready(function() {

        var captureBtn = '%s';
        var masterNm   = '%s';
        var plot1Id    = '%s';
        var plot2Id    = '%s';
        var inputKey   = '%s';

        // Plotly stores the live camera (after user drag/zoom) in the
        // internal scene object, NOT in _fullLayout.scene.camera (which
        // only reflects the initial or last programmatic relayout).
        //
        // The live camera lives at:
        //   el._fullLayout.scene._scene.getCamera()
        // Fallback to _fullLayout.scene.camera if getCamera is unavailable.

        function getLiveCamera(el) {
          // Try the internal WebGL scene first
          try {
            var sceneObj = el._fullLayout.scene._scene;
            if (sceneObj && typeof sceneObj.getCamera === 'function') {
              return sceneObj.getCamera();
            }
          } catch(e) {}

          // Fallback: stale layout camera
          try {
            return el._fullLayout.scene.camera;
          } catch(e) {}

          return null;
        }

        $(document).on('click', '#' + captureBtn, function() {

          var master = $('input[name=\"' + masterNm + '\"]:checked').val();
          var srcId  = (master === 'plot1') ? plot1Id : plot2Id;
          var srcEl  = document.getElementById(srcId);

          if (!srcEl || !srcEl._fullLayout || !srcEl._fullLayout.scene) {
            alert('Master plot not ready. Please wait for it to render.');
            return;
          }

          var cam = getLiveCamera(srcEl);
          if (!cam) {
            alert('Could not read camera from master plot.');
            return;
          }

          var camCopy = JSON.parse(JSON.stringify(cam));

          console.log('Captured LIVE camera:', JSON.stringify(camCopy));

          Shiny.setInputValue(inputKey, camCopy, {priority: 'event'});
        });
      });
    ",
      ns("sc1a2captureBtn"),
      ns("sc1a2master"),
      ns("sc1a2oup1"),
      ns("sc1a2oup2"),
      ns("captured_camera")
    ))),
    
    fluidRow(
      column(
        6, style = "border-right: 2px solid black", h4("Cell information 1"),
        fluidRow(
          column(6,
            selectInput(ns("sc1a2inp1"), "Cell information:",
                        choices = sc1conf$UI, selected = sc1def$meta1)
          ),
          column(6,
            actionButton(ns("sc1a2tog1"), "Customize Plot"),
            conditionalPanel(
              condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2tog1")),
              radioButtons(ns("sc1a2col1"), "Colour (Continuous data):",
                           choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                           selected = "Blue-Yellow-Red"),
              radioButtons(ns("sc1a2ord1"), "Plot order:",
                           choices = c("Max-1st", "Min-1st", "Original", "Random"),
                           selected = "Original", inline = TRUE),
              checkboxInput(ns("sc1a2lab1"), "Show cell info labels", value = TRUE)
            )
          )
        ),
        fluidRow(column(12, uiOutput(ns("sc1a2oup1.ui"))))
      ),
      
      column(
        6, h4("Cell information 2"),
        fluidRow(
          column(6,
            selectInput(ns("sc1a2inp2"), "Cell information:",
                        choices = sc1conf$UI, selected = sc1def$meta2)
          ),
          column(6,
            actionButton(ns("sc1a2tog2"), "Customize plot"),
            conditionalPanel(
              condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a2tog2")),
              radioButtons(ns("sc1a2col2"), "Colour (Continuous data):",
                           choices = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
                           selected = "Blue-Yellow-Red"),
              radioButtons(ns("sc1a2ord2"), "Plot order:",
                           choices = c("Max-1st", "Min-1st", "Original", "Random"),
                           selected = "Original", inline = TRUE),
              checkboxInput(ns("sc1a2lab2"), "Show cell info labels", value = TRUE)
            )
          )
        ),
        fluidRow(column(12, uiOutput(ns("sc1a2oup2.ui"))))
      )
    )
  )
}

########################################## Server ########################################################

scDRcell3D_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    observe_helpers()
    
    # ── 3D reduction choices ──
    dr3d_tbl <- reactive({ scList3DReductions(sc1conf, sc1meta) })
    
    observeEvent(dr3d_tbl(), {
      tbl <- dr3d_tbl()
      choices <- tbl$base
      if (length(choices) == 0) choices <- character(0)
      current <- isolate(input$sc1a2dr3d_base)
      selected <- if (!is.null(current) && current %in% choices) current
                  else if (length(choices) > 0) choices[1] else character(0)
      updateSelectInput(session, "sc1a2dr3d_base", choices = choices, selected = selected)
    }, ignoreInit = FALSE)
    
    # ── Cell filter UI ──
    output$sc1a2sub1.ui <- renderUI({
      sub <- strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(ns("sc1a2sub2"), "Select which cells to show",
                         inline = TRUE, choices = sub, selected = sub)
    })
    
    observeEvent(input$sc1a2sub1non, {
      sub <- strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, "sc1a2sub2", label = "Select which cells to show",
                               choices = sub, selected = NULL, inline = TRUE)
    })
    
    observeEvent(input$sc1a2sub1all, {
      sub <- strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, "sc1a2sub2", label = "Select which cells to show",
                               choices = sub, selected = sub, inline = TRUE)
    })
    
    # ── Plot reactives ──
    plot1 <- reactive({
      scDRcell3D(sc1conf, sc1meta,
        dr3d_base = input$sc1a2dr3d_base, inp1 = input$sc1a2inp1,
        inpsub1 = input$sc1a2sub1, inpsub2 = input$sc1a2sub2,
        inpsiz = input$sc1a2siz, inpcol = input$sc1a2col1,
        inpord = input$sc1a2ord1, inpfsz = input$sc1a2fsz,
        inpasp = input$sc1a2asp, inptxt = input$sc1a2txt, inplab = input$sc1a2lab1)
    })
    
    plot2 <- reactive({
      scDRcell3D(sc1conf, sc1meta,
        dr3d_base = input$sc1a2dr3d_base, inp1 = input$sc1a2inp2,
        inpsub1 = input$sc1a2sub1, inpsub2 = input$sc1a2sub2,
        inpsiz = input$sc1a2siz, inpcol = input$sc1a2col2,
        inpord = input$sc1a2ord2, inpfsz = input$sc1a2fsz,
        inpasp = input$sc1a2asp, inptxt = input$sc1a2txt, inplab = input$sc1a2lab2)
    })
    
    output$sc1a2oup1 <- plotly::renderPlotly({ plot1()$plot })
    output$sc1a2oup1.ui <- renderUI({
      plotly::plotlyOutput(ns("sc1a2oup1"), height = pList[input$sc1a2psz])
    })
    
    output$sc1a2oup2 <- plotly::renderPlotly({ plot2()$plot })
    output$sc1a2oup2.ui <- renderUI({
      plotly::plotlyOutput(ns("sc1a2oup2"), height = pList[input$sc1a2psz])
    })
    
    # ══════════════════════════════════════════════
    # Camera sync: Capture → Display → Apply
    # ══════════════════════════════════════════════
    
    # Store the captured camera (sent from JS via Shiny.setInputValue)
    capturedCam <- reactiveVal(NULL)
    
    observeEvent(input$captured_camera, {
      capturedCam(input$captured_camera)
    })
    
    # Print captured camera values so user can verify
    output$sc1a2camTxt <- renderPrint({
      cam <- capturedCam()
      if (is.null(cam)) {
        cat("No camera captured yet.\nRotate/zoom, then click 'Capture view'.")
      } else {
        cat("Captured camera:\n")
        cat("  eye:    x=", round(cam$eye$x, 3),
                   " y=", round(cam$eye$y, 3),
                   " z=", round(cam$eye$z, 3), "\n")
        cat("  center: x=", round(cam$center$x, 3),
                   " y=", round(cam$center$y, 3),
                   " z=", round(cam$center$z, 3), "\n")
        cat("  up:     x=", round(cam$up$x, 3),
                   " y=", round(cam$up$y, 3),
                   " z=", round(cam$up$z, 3), "\n")
        if (!is.null(cam$projection$type)) {
          cat("  projection:", cam$projection$type, "\n")
        }
      }
    })
    
    # Apply captured camera to the OTHER plot via plotlyProxy
    observeEvent(input$sc1a2applyBtn, {
      cam <- capturedCam()
      if (is.null(cam)) {
        shiny::showNotification("No camera captured yet. Click 'Capture view' first.",
                                type = "warning")
        return()
      }
      
      master   <- input$sc1a2master
      targetId <- if (master == "plot1") "sc1a2oup2" else "sc1a2oup1"
      
      proxy <- plotly::plotlyProxy(targetId, session)
      
      plotly::plotlyProxyInvoke(proxy, "relayout", list(
        scene.camera = list(
          eye    = list(x = cam$eye$x,    y = cam$eye$y,    z = cam$eye$z),
          center = list(x = cam$center$x, y = cam$center$y, z = cam$center$z),
          up     = list(x = cam$up$x,     y = cam$up$y,     z = cam$up$z)
        )
      ))
    })
    
  })
}

############################################### Registration #################################################

register_tab(
  id     = "cellinfo3D_cellinfo3D_multi",
  title  = "CellInfo3D vs CellInfo3D",
  ui     = scDRcell3D_ui,
  server = scDRcell3D_server
)
