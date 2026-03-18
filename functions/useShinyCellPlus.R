useShinyCellPlus <- function(
    shiny.dir, # files from shinycell are
    shinycellplus.dir.src, # modules where shinycellplus 
    rsconnect.deploy = FALSE, # do you want to publish in rsconnect
    data_type = c("RNA", "RNA_ATAC", "SPATIAL"), # what predetermine tabs you want
    enabled_tabs = NULL, # what tabs you want
    overwrite_modules = FALSE, # overwrite modules
    disable_ui_server = TRUE, # this disables the existing ui.R and server.r
    app_title=NULL
) {
  
  shiny.dir <- normalizePath(shiny.dir, mustWork = TRUE)
  shinycellplus.dir.src <- normalizePath(shinycellplus.dir.src, mustWork = TRUE)
  
  message("ShinyCellPlus app generation starting")
  message("Target app directory: ", shiny.dir)
  message("Source ShinyCellPlus directory: ", shinycellplus.dir.src)
  
  # Treat NULL, "", and character(0) as empty
  is_empty <- function(x) {
    is.null(x) || length(x) == 0 || (is.character(x) && all(trimws(x) == ""))
  }
  
  data_type_provided <- !is_empty(data_type)
  tabs_provided <- !is_empty(enabled_tabs)


  
  if (!data_type_provided && !tabs_provided) {
    stop("You must provide either data_type or enabled_tabs. For example data_type can be 'RNA' or RNA_ATAC or SPATIAL. What type of assay do you have?")
  }
  
  
 if (missing(app_title) || is.null(app_title)) {
  stop("App title missing. We are not launching anonymous software today. Please provide a title using app_title='...'.")
 }
  
  
  # Tab catalogue: every   allowed_tabs tab ID per data_type 
  # This is the single source of truth. A tab is only valid for one data_type.
  # Passing a tab that does not belong to the chosen data_type is an error.
  
  module_files <- list.files(file.path(shinycellplus.dir.src, "modules"), recursive = TRUE)

  all_tabs_by_type <- lapply(
    split(module_files, dirname(module_files)),
    function(files) {
      list(
        tab_id   = tools::file_path_sans_ext(basename(files)),
        filename = basename(files)
      )
    }
  )
  
  default_tabs <- NULL
  assays_vec   <- NULL
  
  if (data_type_provided) {
    data_type <- match.arg(data_type, choices = names(all_tabs_by_type))
    
    
    ##?? I am not sure if I am using this for anything
    assays_vec <- switch(
      data_type,
      RNA      = "RNA",
      RNA_ATAC = c("RNA", "ATAC"),
      SPATIAL  = c("Spatial", "RNA")
    )
    ##??
    
    # Default = all   allowed_tabs tabs for this data_type
    default_tabs <- all_tabs_by_type[[data_type]]$tab_id
  }
  
  if(!data_type_provided){ data_type="RNA"
  data_type_provided <- !is_empty(data_type)
  message("You have not provided data_type, we will assume you have Single Cell RNAseq, data_type set to RNA")
  }
  
  if (!tabs_provided) {
    enabled_tabs <- default_tabs
  }
  
  
if (data_type_provided && tabs_provided) {
    # User supplied specific tabs AND a data_type:
    # every requested tab must belong to the   allowed_tabs set for that data_type.
      allowed_tabs      <- all_tabs_by_type[[data_type]]$tab_id
    bad_tabs     <- setdiff(enabled_tabs,   allowed_tabs)

    if (length(bad_tabs) > 0) {
      stop(
        "The following tabs are not valid for data_type = '", data_type, "':\n",
        "  ", paste(bad_tabs, collapse = ", "), "\n\n",
        "  allowed_tabs tabs for '", data_type, "':\n",
        "  ", paste(  allowed_tabs, collapse = ", "), "\n\n",
        "To use tabs from a different data_type, change data_type accordingly.",
        call. = FALSE
      )
    }
    enabled_tabs <-     allowed_tabs 
  } 
  
  message("Enabled tabs : ", paste(enabled_tabs, collapse = ", "))

  if (isTRUE(disable_ui_server)) {
    ui_r <- file.path(shiny.dir, "ui.R")
    server_r <- file.path(shiny.dir, "server.R")
    
    if (file.exists(ui_r) || file.exists(server_r)) {
      warning(
        paste(
          "ui.R and or server.R detected in the app directory.",
          "Shiny will prioritise these files over app.R.",
          "To ensure the modular ShinyCellPlus app is used,",
          "ui.R and server.R will be disabled by renaming them.",
          "Backup files with extension .bak will be created."
        ),
        call. = FALSE
      )
    }
    
    if (file.exists(ui_r)) {
      file.rename(ui_r, file.path(shiny.dir, "ui.R.bak"))
      message("Renamed ui.R to ui.R.bak")
    }
    
    if (file.exists(server_r)) {
      file.rename(server_r, file.path(shiny.dir, "server.R.bak"))
      message("Renamed server.R to server.R.bak")
    }
  }
  
  
  ###### Here I am trying to copy scripts directly by name AQUI VOYYY
  src_modules <- file.path(shinycellplus.dir.src, "modules/",data_type,all_tabs_by_type[[data_type]]$filename)
  src_modules_dir <- file.path(shinycellplus.dir.src, "modules/")
  dst_modules <- file.path(shiny.dir, "modules/")
  
  if (!dir.exists(src_modules_dir)) {
    stop("Could not find 'modules' folder in shinycellplus.dir.src: ", src_modules_dir)
  }
  
  if (dir.exists(dst_modules) && isTRUE(overwrite_modules)) {
    warning(
      paste(
        "An existing modules directory was found.",
        "overwrite_modules = TRUE, so it will be removed and replaced.",
        "Any local modifications inside modules/ will be lost."
      ),
      call. = FALSE
    )
  }
  
  if (!dir.exists(dst_modules) || isTRUE(overwrite_modules)) {
    if (dir.exists(dst_modules)) unlink(dst_modules, recursive = TRUE, force = TRUE)
    dir.create(dst_modules, recursive = TRUE, showWarnings = FALSE)
   
     ok <- file.copy(src_modules, dst_modules, recursive = FALSE)
    
    failed <- src_modules[!ok]
    if (length(failed) > 0) {
      stop("Failed to copy the following files to ", dst_modules, ":\n  ",
           paste(basename(failed), collapse = "\n  "))
    }
    
    message("Copied ", sum(ok), " module(s) into: ", dst_modules)
  } else {
    message("Using existing modules/ folder in: ", dst_modules)
  }
  
  dir_inputs <- shiny.dir
  
  assays_str <- paste0("c(", paste(sprintf('"%s"', assays_vec), collapse = ", "), ")")
  enabled_tabs_str <- paste0("c(", paste(sprintf('"%s"', enabled_tabs), collapse = ", "), ")")
  
  app_modules <- sprintf(
    '## Auto generated by useShinyCellPlus

message("Starting ShinyCellPlus modular app")

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
library(shinythemes)
library(shinydashboard)
library(tidyverse)
library(sortable)
library(plotly)
library(FlexDotPlot) #devtools::install_github("Simon-Leonard/FlexDotPlot")
library(RColorBrewer)
library(ggforce)
library(limma) #BiocManager::install("limma")
library(edgeR) #BiocManager::install("edgeR")



### Useful stuff 
# Colour palette 
cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154")) 
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple") 
 
# Panel sizes 
pList = c("400px", "600px", "800px") 
names(pList) = c("Small", "Medium", "Large") 
pList2 = c("500px", "700px", "900px") 
names(pList2) = c("Small", "Medium", "Large") 
pList3 = c("600px", "800px", "1000px") 
names(pList3) = c("Small", "Medium", "Large") 
sList = c(18,24,30) 
names(sList) = c("Small", "Medium", "Large") 
lList = c(5,6,7) 
names(lList) = c("Small", "Medium", "Large") 
 
# Function to extract legend 
g_legend <- function(a.gplot){  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
  legend <- tmp$grobs[[leg]]  
  legend 
}  
 
# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
  oupTheme = theme( 
    text =             element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line =   element_line(colour = "black"), 
    axis.ticks =  element_line(colour = "black", size = base_size / 20), 
    axis.title =  element_text(face = "bold"), 
    axis.text =   element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key =      element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 
 
app_title <- "%1$s"

dir_inputs <- "%2$s/"

sc1conf <- readRDS(file.path(dir_inputs, "sc1conf.rds"))
sc1def  <- readRDS(file.path(dir_inputs, "sc1def.rds"))
sc1gene <- readRDS(file.path(dir_inputs, "sc1gene.rds"))
sc1meta <- readRDS(file.path(dir_inputs, "sc1meta.rds"))

if (file.exists(file.path(dir_inputs, "markergenes_lists.parquet"))) {
  markers_list <- file.path(dir_inputs, "markergenes_lists.parquet")
} else {
  markers_list <- NULL
}

assays <- %3$s
assays_vec <- unique(sc1conf$assay)

tab_registry <- list()

register_tab <- function(id, title, ui, server) {

  fn_args <- names(formals(server))
  if (is.null(fn_args)) fn_args <- character()

  has_id <- "id" %%in%% fn_args
  has_input <- "input" %%in%% fn_args
  has_output <- "output" %%in%% fn_args
  has_session <- "session" %%in%% fn_args

  if (!isTRUE(has_id)) {
    stop("register_tab: server function must have argument id for tab: ", id)
  }

  if (isTRUE(has_input) || isTRUE(has_output) || isTRUE(has_session)) {
    warning(
      paste(
        "Tab", id, "server declares input and or output and or session as arguments.",
        "For Shiny module style, server should be function(id, ...) and call moduleServer inside.",
        "This tab will still be registered, but consider updating the signature."
      ),
      call. = FALSE
    )
  }

  tab_registry[[id]] <<- list(
    title  = title,
    ui     = ui,
    server = server
  )
}

get_tab_ids <- function(enabled_tabs = NULL) {
  all_ids <- names(tab_registry)
  if (is.null(enabled_tabs) || length(enabled_tabs) == 0) all_ids else intersect(enabled_tabs, all_ids)
}

get_app_dir <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && is.character(ofile) && nzchar(ofile)) {
    return(normalizePath(dirname(ofile), winslash = "/", mustWork = TRUE))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    f <- sub("^--file=", "", file_arg)
    if (nzchar(f) && file.exists(f)) {
      return(normalizePath(dirname(f), winslash = "/", mustWork = TRUE))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

app_dir <- get_app_dir()
modules_dir <- file.path(app_dir, "modules")

if (!dir.exists(modules_dir)) {
  stop("Modules dir not found: ", modules_dir, " Current working directory: ", getwd())
}

for (f in list.files(modules_dir, full.names = TRUE, pattern = "\\\\.[Rr]$")) {
  message("Sourcing module: ", basename(f))
  tryCatch(
    source(f, local = environment()),
    error = function(e) {
      warning(
        paste0("Skipping module due to error: ", basename(f), " | ", conditionMessage(e)),
        call. = FALSE
      )
    }
  )
}


if (length(tab_registry) == 0) {
  warning(
    paste(
      "No modules registered any tabs.",
      "This usually means module files were not sourced correctly",
      "or register_tab() was not called.",
      "Check modules/ and module file names."
    ),
    call. = FALSE
  )
}

enabled_tabs <-  %4$s

missing_tabs <- setdiff(enabled_tabs, names(tab_registry))
if (length(missing_tabs) > 0) {
  warning(
    paste(
      "The following requested tabs were not found:",
      paste(missing_tabs, collapse = ", "),
      "They will be ignored."
    ),
    call. = FALSE
  )
}

tab_ids <- get_tab_ids(enabled_tabs)

tab_panels <- lapply(tab_ids, function(k) {
  tabPanel(
    tab_registry[[k]]$title,
    tab_registry[[k]]$ui(id = k, sc1conf = sc1conf, sc1def = sc1def)
  )
})

ui <- fluidPage( theme = shinytheme("cerulean"),
      tags$head(
        tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}")),
        tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))
      ),
      titlePanel(app_title),
      do.call(navbarPage, c(list(title = NULL), tab_panels)),
      tags$hr(),
      tags$p(
        style = "font-size: 90%%; color: #666;",
        em(
          "This application was generated using ShinyCellPlus. ",
          "Tabs are dynamically loaded from modular components.",
          "Monash Genomics and Bioinformatics Platform. ShinyCellPlus: ShinyCell Package Customized by MGBP v.1 Date: Jan 2026"
        )
      ),
      br(), br(), br(), br(), br()
    )
  


server <- function(input, output, session) {
  lapply(tab_ids, function(k) {

    srv <- tab_registry[[k]]$server

    args_to_pass <- list(
      id = k,
      sc1conf = sc1conf,
      sc1meta = sc1meta,
      sc1gene = sc1gene,
      sc1def  = sc1def,
      markers_list = markers_list,
      assays = assays,
      dir_inputs = dir_inputs
    )

    keep <- intersect(names(args_to_pass), names(formals(srv)))
    do.call(srv, args_to_pass[keep])
  })
}


shinyApp(ui, server)
',app_title,
dir_inputs,
assays_str,
enabled_tabs_str
  )



app_path <- file.path(shiny.dir, "app.R")
writeLines(app_modules, con = app_path)
message("Wrote app.R to: ", app_path)

if (isTRUE(rsconnect.deploy)) {
  library(rsconnect)
  library(jsonlite)
  
  rsconnect::writeManifest(appDir = shiny.dir)
  message("Wrote rsconnect manifest in: ", shiny.dir)
  
  dir_prefix <- shiny.dir
  manifest_path <- file.path(shiny.dir, "manifest.json")
  m <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  
  stopifnot(!is.null(m$files))
  old_keys <- names(m$files)
  
  is_target <- grepl("\\.(rds|h5||parquet)$", old_keys, ignore.case = TRUE)
  
  targets <- old_keys[is_target]
  if (length(targets) == 0) stop("No .rds or .h5 files found in manifest$files")
  
  for (k in targets) {
    new_k <- paste0(dir_prefix,"/", k)
    if (!is.null(m$files[[new_k]])) stop("Target key already exists: ", new_k)
    m$files[[new_k]] <- m$files[[k]]
    m$files[[k]] <- NULL
  }
  
  writeLines(toJSON(m, pretty = TRUE, auto_unbox = TRUE), manifest_path)
  
  cat("Updated keys:\n")
  cat(paste0("  ", targets, " -> ", paste0(dir_prefix, targets)), sep = "\n")
  cat("\n")
  
}

invisible(app_path)
}
