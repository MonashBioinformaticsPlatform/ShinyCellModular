# Developer Instructions — Creating New Modules

## Overview

Each tab in ShinyCellModular is a self-contained R file called a **module**. Modules live inside `modules/<data_type>/` and are sourced automatically by the app at startup. No changes to `app.R` or `useShinyCellModular()` are needed when adding a new module — dropping a file into the right subfolder is enough.

Each module file contains exactly four sections in order:

1. **Functions** — pure R functions that do the computation and plotting
2. **UI function** — defines the tab layout and inputs
3. **Server function** — defines the reactive logic, renders plots and tables
4. **Registration** — registers the tab with the app via `register_tab()`

---

## Folder structure

```
modules/
  RNA/
    bubble_heatmap.R
    cellinfo_cellinfo.R
    ...
  RNA_ATAC/
    bubble_heatmap_multi.R
    coverage_plot.R
    ...
  CROPSEQ/
  SPATIAL/
develop/
  my_new_tab_in_progress.R
```

- `modules/` — production-ready tabs, loaded automatically
- `develop/` — work in progress, not loaded by the app
- Each `data_type` subfolder corresponds to the `data_type` argument in `useShinyCellModular()`

When you are ready to test a module in the app, move it from `develop/` into the appropriate `modules/<data_type>/` subfolder.

---

## Naming conventions

| What | Convention | Example |
|---|---|---|
| Module filename | `tab_id.R` | `my_analysis.R` |
| Tab id | lowercase, underscores | `my_analysis` |
| UI function | `<descriptive_name>_ui` | `myAnalysis_ui` |
| Server function | `<descriptive_name>_server` | `myAnalysis_server` |
| Input/output ids | `sc1<letter><number><name>` | `sc1f1grp`, `sc1f1oup` |
| RNA_ATAC variants | `<tab_id>_multi.R` | `my_analysis_multi.R` |

The `sc1<letter><number>` prefix system avoids input id clashes across modules. Check existing modules to pick a letter/number combination that is not already in use. The `createSCModuleTemplate()` helper does this automatically.

---

## Quick start — using the template helper

The safest way to start a new module is with the helper function, which scans existing modules and avoids prefix and id clashes:

```r
source("functions/createSCModuleTemplate.R")

createSCModuleTemplate(
  module_dir = "ShinyCellModular/modules/RNA/",
  tab_id     = "my_analysis",
  tab_title  = "My Analysis"
)
```

This creates a scaffolded file with the correct structure. Move it to `develop/` if you want to work on it before it is live.

---

## Section 1 — Functions

Place all computation, data wrangling, and plotting logic here as plain R functions. Keep them independent of Shiny — they should work outside the app and be testable on their own.

```r
############################################### Functions ###########################################

# Main plot function — takes pre-processed inputs, returns a ggplot or grob
myAnalysis_plot <- function(inpConf, inpMeta, inpGrp, inpGene,
                             inpcols, inpfsz, save = FALSE) {

  # ... computation ...

  ggOut <- ggplot2::ggplot(...) + ...

  ggLeg <- g_legend(ggOut)
  ggOut <- ggOut + ggplot2::theme(legend.position = "none")

  if (!isTRUE(save)) {
    gridExtra::grid.arrange(ggOut, ggLeg, heights = c(7, 2),
                            layout_matrix = rbind(c(1), c(2)))
  } else {
    gridExtra::arrangeGrob(ggOut, ggLeg, heights = c(7, 2),
                           layout_matrix = rbind(c(1), c(2)))
  }
}
```

**Tips:**
- Pass `save = FALSE` for on-screen rendering, `save = TRUE` for `ggsave()` — `arrangeGrob` does not draw to screen, `grid.arrange` does
- Use `shiny::validate(shiny::need(...))` inside plot functions to show user-friendly error messages instead of crashing
- Use `g_legend()` (available globally in `app.R`) to extract and reattach the legend separately, which gives you control over layout
- Use `sctheme()` (available globally) for consistent plot styling
- Use `cList[[inpcols]]` and `sList[inpfsz]` for colour scheme and font size, which match the app-wide options

---

## Section 2 — UI function

The UI function defines everything the user sees in the tab: inputs on the left, plot and download buttons on the right.

```r
################################################# UI #################################################

myAnalysis_ui <- function(id, sc1conf, sc1def) {

  ns <- shiny::NS(id)                         # always the first line

  tabPanel(
    HTML("My Analysis"),                      # tab label shown in the navbar
    h4("My analysis title"),
    "Brief description of what this tab does.", br(), br(),

    fluidRow(
      column(
        3, style = "border-right: 2px solid black",

        # ── Inputs ──
        selectInput(
          ns("sc1f1grp"), "Group by:",
          choices  = sc1conf[grp == TRUE]$UI,
          selected = sc1conf[grp == TRUE]$UI[1]
        ),

        # ── Toggle sections ──
        actionButton(ns("sc1f1tog"), "Customize plot"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1f1tog")),
          radioButtons(
            ns("sc1f1cols"), "Colour scheme:",
            choices  = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple"),
            selected = "Blue-Yellow-Red"
          ),
          radioButtons(
            ns("sc1f1psz"), "Plot size:",
            choices  = c("Small", "Medium", "Large"),
            selected = "Medium", inline = TRUE
          ),
          radioButtons(
            ns("sc1f1fsz"), "Font size:",
            choices  = c("Small", "Medium", "Large"),
            selected = "Medium", inline = TRUE
          )
        )
      ),

      column(
        9,
        h4(htmlOutput(ns("sc1f1oupTxt"))),
        uiOutput(ns("sc1f1oup.ui")),
        downloadButton(ns("sc1f1oup.pdf"), "Download PDF"),
        downloadButton(ns("sc1f1oup.png"), "Download PNG"),
        br(),
        div(
          style = "display:inline-block",
          numericInput(ns("sc1f1oup.h"), "PDF / PNG height:", width = "138px",
                       min = 4, max = 20, value = 10, step = 0.5)
        ),
        div(
          style = "display:inline-block",
          numericInput(ns("sc1f1oup.w"), "PDF / PNG width:", width = "138px",
                       min = 4, max = 20, value = 10, step = 0.5)
        )
      )
    )
  )
}
```

### UI rules and gotchas

**Always wrap input ids with `ns()`**
Every `inputId`, `outputId`, and element `id` inside the UI function must be wrapped in `ns()`. Missing this is the most common source of broken modules.

```r
# correct
selectInput(ns("sc1f1grp"), ...)
plotOutput(ns("sc1f1oup"))

# wrong — will conflict with inputs from other modules
selectInput("sc1f1grp", ...)
```

**`conditionalPanel` conditions use the namespaced id inside a sprintf**
The condition string is evaluated as JavaScript, and the input name must include the namespace prefix. Always use this pattern:

```r
conditionalPanel(
  condition = sprintf("input['%s'] %% 2 == 1", ns("sc1f1tog")),
  ...
)
```

Note the `%%` to escape the modulo operator inside `sprintf`.

**Use `uiOutput` + `renderUI` for the plot panel**
This allows the plot height to be reactive (driven by the size selector):

```r
# in UI
uiOutput(ns("sc1f1oup.ui"))

# in server
output$sc1f1oup.ui <- renderUI({
  plotOutput(ns("sc1f1oup"), height = pList3[input$sc1f1psz])
})
```

**Available global size lists**

| Object | Small | Medium | Large | Use for |
|---|---|---|---|---|
| `pList` | 400px | 600px | 800px | standard plots |
| `pList2` | 500px | 700px | 900px | slightly taller plots |
| `pList3` | 600px | 800px | 1000px | heatmaps / tall gene lists |
| `sList` | 18 | 24 | 30 | font size (numeric, passed to `sctheme`) |

**Selectize gene inputs use `selectizeInput` not `selectInput`**
Gene name inputs must use server-side selectize for performance with large gene lists (see server setup below).

---

## Section 3 — Server function

The server function handles all reactive logic: updating inputs, rendering plots, and wiring download handlers.

```r
############################################## Server ################################################

myAnalysis_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  shiny::moduleServer(id, function(input, output, session) {

    ns <- session$ns                          # always the first line

    # ── Standard setup (keep at the top of every server) ──
    observe_helpers()
    optCrt <- "{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"

    # Server-side selectize for gene inputs (only include the ones your tab uses)
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
                         choices = c(sc1conf[is.na(fID)]$UI, names(sc1gene)),
                         selected = sc1conf[is.na(fID)]$UI[1],
                         options = list(
                           maxOptions = length(sc1conf[is.na(fID)]$UI) + 3,
                           create = TRUE, persist = TRUE, render = I(optCrt)))

    # ── Output text ──
    output$sc1f1oupTxt <- renderUI({
      HTML("Some status text or summary shown above the plot")
    })

    # ── Render plot ──
    output$sc1f1oup <- renderPlot({
      p <- myAnalysis_plot(
        sc1conf, sc1meta,
        input$sc1f1grp,
        input$sc1f1cols, input$sc1f1fsz
      )
      if (inherits(p, "grob") || inherits(p, "gtable")) {
        grid::grid.newpage()
        grid::grid.draw(p)
      } else {
        print(p)
      }
    })

    # ── Plot panel size ──
    output$sc1f1oup.ui <- renderUI({
      plotOutput(ns("sc1f1oup"), height = pList3[input$sc1f1psz])
    })

    # ── Download PDF ──
    output$sc1f1oup.pdf <- downloadHandler(
      filename = function() {
        paste0("sc1myAnalysis_", input$sc1f1grp, ".pdf")
      },
      content = function(file) {
        ggplot2::ggsave(
          file, device = "pdf",
          height = input$sc1f1oup.h, width = input$sc1f1oup.w,
          plot = myAnalysis_plot(
            sc1conf, sc1meta,
            input$sc1f1grp,
            input$sc1f1cols, input$sc1f1fsz,
            save = TRUE
          )
        )
      }
    )

    # ── Download PNG ──
    output$sc1f1oup.png <- downloadHandler(
      filename = function() {
        paste0("sc1myAnalysis_", input$sc1f1grp, ".png")
      },
      content = function(file) {
        ggplot2::ggsave(
          file, device = "png",
          height = input$sc1f1oup.h, width = input$sc1f1oup.w,
          plot = myAnalysis_plot(
            sc1conf, sc1meta,
            input$sc1f1grp,
            input$sc1f1cols, input$sc1f1fsz,
            save = TRUE
          )
        )
      }
    )

  })
}
```

### Server rules and gotchas

**Always initialise `ns <- session$ns` at the top**
This is needed wherever you generate UI elements inside the server (e.g. inside `renderUI()`). Without it, output ids will not be namespaced correctly and the plot will not render.

```r
# correct — use ns() when declaring UI elements inside renderUI
output$sc1f1oup.ui <- renderUI({
  plotOutput(ns("sc1f1oup"), height = pList3[input$sc1f1psz])
})

# also correct — reading inputs does not need ns()
input$sc1f1grp
```

**Server-side selectize setup**
Only include `updateSelectizeInput()` calls for the gene inputs your module actually uses. The block is boilerplate from the original ShinyCell structure — it is safe to remove entries your module does not need. Keep `observe_helpers()` regardless.

**Reading from HDF5**
Gene expression data is stored in `sc1gexpr.h5` inside `dir_inputs`. Access it like this:

```r
h5file <- hdf5r::H5File$new(file.path(dir_inputs, "sc1gexpr.h5"), mode = "r")
on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
h5data <- h5file[["grp"]][["data"]]
val <- h5data$read(args = list(sc1gene[iGene], quote(expr=)))
```

Always use `on.exit` to close the file — if the reactive errors before reaching `h5file$close_all()`, the file handle will otherwise remain open.

**Reading raw counts (pseudobulk tabs)**
Raw counts are in `sc1counts.h5` inside `dir_inputs`. This file is only present if `do_counts_h5 = TRUE` was used in `prepShinyCellModular()`. Check for its existence before using it:

```r
counts_path <- file.path(dir_inputs, "sc1counts.h5")
shiny::validate(shiny::need(file.exists(counts_path),
  "sc1counts.h5 not found. Re-run prepShinyCellModular with do_counts_h5 = TRUE."))
```

**Rendering grobs**
When your plot function returns a `gridExtra` layout (a grob or gtable), use `grid.draw()` rather than `print()`:

```r
if (inherits(p, "grob") || inherits(p, "gtable")) {
  grid::grid.newpage()
  grid::grid.draw(p)
} else {
  print(p)
}
```

**Available global objects inside the server**
These are passed in by `useShinyCellModular()` and available as arguments if declared in the server signature:

| Object | Type | Contents |
|---|---|---|
| `sc1conf` | `data.table` | Metadata configuration — column IDs, UI labels, factor levels |
| `sc1meta` | `data.table` | Per-cell metadata — one row per cell |
| `sc1gene` | named vector | Gene name → HDF5 index mapping |
| `sc1def` | list | Default values — `gene1`, `gene2`, `grp1`, `grp2`, `genes` |
| `dir_inputs` | string | Path to the prepared data folder |
| `markers_list` | string or NULL | Path to `markergenes_lists.parquet`, or NULL if absent |
| `sc1conf_atac` | `data.table` or NULL | ATAC metadata config, NULL if no ATAC data |
| `sc1meta_atac` | `data.table` or NULL | ATAC per-cell metadata |
| `sc1gene_atac` | named vector or NULL | ATAC gene/peak index mapping |
| `sc1fragmentpaths` | list or NULL | Fragment file paths and barcodes |
| `sc1annotation` | GRanges or NULL | Genome annotation |
| `sc1peaks` | GRanges or NULL | Peak ranges |
| `sc1links` | GRanges or NULL | Peak-to-gene links |

Only declare the arguments your module needs — the app wires them via `formals()` matching, so undeclared arguments are simply not passed.

**Accessing metadata columns**
Use `sc1conf` to map UI labels to internal column IDs:

```r
# get the data.table column name for a UI label
col_id <- sc1conf[UI == input$sc1f1grp]$ID

# get the factor levels for a grouping variable
levels_str <- sc1conf[UI == input$sc1f1grp]$fID   # pipe-separated string
levels_vec <- strsplit(levels_str, "\\|")[[1]]

# access the actual values
sc1meta[[col_id]]
```

---

## Section 4 — Registration

Every module must end with a `register_tab()` call. This is what makes the tab discoverable by the app.

```r
############################################### Registration ###########################################

register_tab(
  id     = "my_analysis",     # must match the filename without .R extension
  title  = "My Analysis",     # label shown in the navbar
  ui     = myAnalysis_ui,
  server = myAnalysis_server
)
```

The `id` must:
- Match the filename without the `.R` extension exactly
- Be unique across all modules in all subfolders
- Use lowercase and underscores only

---

## Dependencies

If your module uses additional packages, add them in two places:

1. **`renv`** — run `renv::snapshot()` after installing
2. **`useShinyCellModular()`** — add the `library()` call to the `template_app` string inside `useShinyCellModular.R` so it is included in the generated `app.R`

Current dependencies available globally in the app:

```r
library(shiny)
library(shinyhelper)
library(shinyjs)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(hdf5r)
library(ggdendro)
library(gridExtra)
library(arrow)
library(rsconnect)
library(shinythemes)
library(shinydashboard)
library(tidyverse)
library(sortable)
library(plotly)
library(FlexDotPlot)
library(RColorBrewer)
library(ggforce)
library(limma)
library(edgeR)
library(gtable)
```

---

## Checklist for a new module

- [ ] Filename matches the tab id (`my_analysis.R` → `id = "my_analysis"`)
- [ ] File is in the correct `modules/<data_type>/` subfolder
- [ ] All input and output ids are wrapped in `ns()` in the UI
- [ ] `conditionalPanel` conditions use `sprintf("input['%s'] ...", ns(...))`
- [ ] Server starts with `ns <- session$ns`
- [ ] `observe_helpers()` is called at the top of the server
- [ ] HDF5 files are closed with `on.exit(h5file$close_all())`
- [ ] Plot function accepts `save = FALSE` and returns either a grob (via `arrangeGrob`) or a ggplot
- [ ] `renderPlot` handles both grob and ggplot return types
- [ ] `register_tab()` id matches the filename exactly
- [ ] Any new package dependencies are added to `renv` and `useShinyCellModular.R`
- [ ] Tested in `develop/` before moving to `modules/`
