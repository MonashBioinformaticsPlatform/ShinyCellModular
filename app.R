########################################### Load libraries ###########################################

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

########################################### Load files ###########################################

dir_inputs <- "<FULL_PATH_WILL_BE_INJECTED_HERE>/"

sc1conf <- readRDS(file.path(dir_inputs, "sc1conf.rds"))
sc1def  <- readRDS(file.path(dir_inputs, "sc1def.rds"))
sc1gene <- readRDS(file.path(dir_inputs, "sc1gene.rds"))
sc1meta <- readRDS(file.path(dir_inputs, "sc1meta.rds"))

markers_list      <- read_parquet(file.path(dir_inputs, "markergenes_lists.parquet"))

########################################### Style and helpers ###########################################

cList <- list(
  c("grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
    "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"),
  c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
    "#FEE090", "#FDAE61", "#F46D43", "#D73027")[c(1, 1:9, 9)],
  c("#FDE725", "#AADC32", "#5DC863", "#27AD81", "#21908C",
    "#2C728E", "#3B528B", "#472D7B", "#440154")
)
names(cList) <- c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")

pList  <- c("400px", "600px", "800px");  names(pList)  <- c("Small", "Medium", "Large")
pList2 <- c("500px", "700px", "900px");  names(pList2) <- c("Small", "Medium", "Large")
pList3 <- c("600px", "800px", "1000px"); names(pList3) <- c("Small", "Medium", "Large")
sList  <- c(18, 24, 30);                 names(sList)  <- c("Small", "Medium", "Large")
lList  <- c(5, 6, 7);                    names(lList)  <- c("Small", "Medium", "Large")

g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5) {
  oupTheme <- theme(
    text            = element_text(size = base_size, family = "Helvetica"),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line       = element_line(colour = "black"),
    axis.ticks      = element_line(colour = "black", size = base_size / 20),
    axis.title      = element_text(face = "bold"),
    axis.text       = element_text(size = base_size),
    axis.text.x     = element_text(angle = Xang, hjust = XjusH),
    legend.position = "bottom",
    legend.key      = element_rect(colour = NA, fill = NA)
  )
  if (!XYval) {
    oupTheme <- oupTheme + theme(
      axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y  = element_blank(), axis.ticks.y = element_blank()
    )
  }
  oupTheme
}



########################################### Registry in memory ###########################################

tab_registry <- list()

register_tab <- function(id, title, ui, server) {
  tab_registry[[id]] <<- list(
    title  = title,
    ui     = ui,
    server = server
  )
}

get_tab_ids <- function(enabled_tabs = NULL) {
  all_ids <- names(tab_registry)
  if (is.null(enabled_tabs) || length(enabled_tabs) == 0) {
    all_ids
  } else {
    intersect(enabled_tabs, all_ids)
  }
}


################################################ Modules ############################################################

## modules directory lives inside this app folder
modules_dir <- file.path(getwd(), "code", "modules")
for (f in list.files(modules_dir, full.names = TRUE)) {
  source(f, local = TRUE)
}

########################################### Choose which tabs to show ###################################

## default: all registered modules
enabled_tabs <- NULL   ## or eg c("cellinfo_geneexpr", "umap", "viobox")

tab_ids <- get_tab_ids(enabled_tabs)


enabled_tabs <- c("umap", "viobox")  # keys in tab_registry

########################################### UI ###########################################

ui <- navbarPage(
  "ShinyCell modular app",
  tabPanel(
    title = "App",
    fluidPage(
      shinyjs::useShinyjs(),
      tags$head(
        tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}")),
        tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))
      ),
      titlePanel("ShinyCell Customized MGBP"),

      # inner navbar with modular tabs
      navbarPage(
        NULL,
        !!!lapply(enabled_tabs, function(k) {
          tabPanel(
            tab_registry[[k]]$title,
            tab_registry[[k]]$ui(id = k, sc1conf = sc1conf, sc1def = sc1def)
          )
        })
      ),

      # footer
      br(),
      p("", style = "font-size: 125%;"),

      p(
        em("This webpage was made using a customized version of:"),
        br(),
        tags$ul(
          tags$li(
            "ShinyCell ",
            a("GitHub",
              href   = "https://github.com/SGDDNB/ShinyCell",
              target = "_blank"
            )
          ),
          tags$li(
            "ShinyCellPlus (MGBP) ",
            a("GitHub",
              href   = "https://github.com/MonashBioinformaticsPlatform/ShinycellPlus",
              target = "_blank"
            )
          )
        )
      ),

      br(), br(), br(), br(), br()
    )
  )
)


########################################### Server ###########################################

server <- function(input, output, session) {
  lapply(enabled_tabs, function(k) {
    tab_registry[[k]]$server(
      id      = k,
      sc1conf = sc1conf,
      sc1meta = sc1meta,
      sc1gene = sc1gene,
      sc1def  = sc1def,
      markers_list     = markers_list,
      markers_list_roc = markers_list_roc,
      assays           = assays
    )
  })
}


########################################### Run app ###########################################
shinyApp(ui, server)





