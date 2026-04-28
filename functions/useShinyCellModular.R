useShinyCellModular <- function(
    shiny.dir, # files from shinycell are
    shinycellmodular.dir.src, # modules where shinycellmodular 
    rsconnect.deploy = FALSE, # do you want to publish in rsconnect
    data_type = c("RNA", "RNA_ATAC", "SPATIAL"), # what predetermine tabs you want
    enabled_tabs = NULL, # what tabs you want
    overwrite_modules = FALSE, # overwrite modules
    disable_ui_server = TRUE, # this disables the existing ui.R and server.r
    app_title=NULL
) {
  
  shiny.dir <- normalizePath(shiny.dir, mustWork = TRUE)
  shinycellmodular.dir.src <- normalizePath(shinycellmodular.dir.src, mustWork = TRUE)
  
  message("ShinyCellModular app generation starting")
  message("Target app directory: ", shiny.dir)
  message("Source ShinyCellModular directory: ", shinycellmodular.dir.src)
  
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
  
  module_files <- list.files(file.path(shinycellmodular.dir.src, "modules"), recursive = TRUE)

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
      selected_tabs<- intersect(allowed_tabs,enabled_tabs)

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
    enabled_tabs <-  selected_tabs 
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
          "To ensure the modular ShinyCellModular app is used,",
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
  
  
  idx_tabs<-which(all_tabs_by_type[[data_type]]$tab_id %in% enabled_tabs)
  
  # copy exact requested tabs only, avoid copying excessive amount of files if this app grow
  src_modules <- file.path(shinycellmodular.dir.src, "modules/",data_type,all_tabs_by_type[[data_type]]$filename[idx_tabs])
  src_modules_dir <- file.path(shinycellmodular.dir.src, "modules/")
  dst_modules <- file.path(shiny.dir, "modules/")
  
  if (!dir.exists(src_modules_dir)) {
    stop("Could not find 'modules' folder in shinycellmodular.dir.src: ", src_modules_dir)
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
  
  assays_str <- paste0("c(", paste(sprintf('"%s"', data_type), collapse = ", "), ")")
 # assays_str <- data_type
  enabled_tabs_str <- paste0("c(", paste(sprintf('"%s"', enabled_tabs), collapse = ", "), ")")
  
template_app<-  '## Auto generated by useShinyCellModular

message("Starting ShinyCellModular modular app")

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
library(ggseqlogo)


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
  leg <- which(sapply(tmp$grobs, function(x){x$name}) == "guide-box")  
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
 
app_title <- "__APP_TITLE__"

dir_inputs <- "__DIR_INPUTS__/"

sc1conf <- readRDS(file.path(dir_inputs, "sc1conf.rds"))
sc1def  <- readRDS(file.path(dir_inputs, "sc1def.rds"))
sc1gene <- readRDS(file.path(dir_inputs, "sc1gene.rds"))
sc1meta <- readRDS(file.path(dir_inputs, "sc1meta.rds"))

# there must be a better way to do this for all alternative assays

if (file.exists(file.path(dir_inputs,"ATAC"))) {
  atac_dir         <- file.path(dir_inputs, "ATAC")
  sc1conf_atac     <- tryCatch({readRDS(file.path(atac_dir, "sc1conf.rds"))},error = function(e) {return(NULL)})
  sc1def_atac      <- tryCatch({readRDS(file.path(atac_dir, "sc1def.rds"))},error = function(e){return(NULL)})
  sc1gene_atac     <- tryCatch({readRDS(file.path(atac_dir, "sc1gene.rds"))},error = function(e){return(NULL)})
  sc1meta_atac     <- tryCatch({readRDS(file.path(atac_dir, "sc1meta.rds"))},error = function(e){return(NULL)})
  sc1fragmentpaths <- tryCatch({readRDS(file.path(atac_dir, "sc1fragmentpaths.rds"))},error = function(e){return(NULL)})
  sc1annotation    <- tryCatch({readRDS(file.path(atac_dir, "sc1annotation.rds"))},error = function(e){return(NULL)})
  sc1peaks         <- tryCatch({readRDS(file.path(atac_dir, "sc1peaks.rds"))},error = function(e){return(NULL)})
  sc1links         <- tryCatch({readRDS(file.path(atac_dir, "sc1links.rds"))},error = function(e) {return(NULL)})
} else { sc1conf_atac     <- NULL
          sc1def_atac      <- NULL
          sc1gene_atac     <- NULL
          sc1meta_atac     <- NULL
          sc1fragmentpaths <- NULL
          sc1annotation    <- NULL
          sc1peaks         <- NULL
          sc1links         <- NULL
        }

if (file.exists(file.path(dir_inputs, "markergenes_lists.parquet"))) {
  markers_list <- file.path(dir_inputs, "markergenes_lists.parquet")
} else {
  markers_list <- NULL 
}

assays <- __ASSAYS__ # still unclear if I am using this for anything
assays_vec <- unique(sc1conf$assay)

tab_registry <- list()

register_tab <- function(id, title, ui, server) {

  fn_args <- names(formals(server))
  if (is.null(fn_args)) fn_args <- character()

  has_id <- "id" %in% fn_args
  has_input <- "input" %in% fn_args
  has_output <- "output" %in% fn_args
  has_session <- "session" %in% fn_args

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
  ofile <- tryCatch({sys.frame(1)$ofile}, error = function(e) {retunr(NULL)})
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
    {source(f, local = environment())},
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

enabled_tabs <- __ENABLED_TABS__

missing_tabs <- setdiff(enabled_tabs, names(tab_registry))
if (length(missing_tabs) > 0) {
  warning(
    paste(
      "The following requested tabs were not found:",
      paste(missing_tabs, collapse = ", "),
      "They will be ignored. Make sure the registry in the modules matches your requested tab"
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
        tags$style(HTML(".navbar{min-height:36px;} .navbar-default .navbar-nav>li>a{padding-top:8px;padding-bottom:8px;font-size:13px;font-weight:bold;} .navbar-default .navbar-brand{padding-top:8px;padding-bottom:8px;font-size:13px;font-weight:bold;border-right:1px solid rgba(255,255,255,0.3);margin-right:4px;} .navbar-collapse{padding-top:0;padding-bottom:0;}"))
      ),
      do.call(navbarPage, c(list(title = app_title), tab_panels)),
      tags$hr(),
      tags$p(
        style = "font-size: 90%; color: #666;",
        em(
          "This application was generated using ShinyCellModular. ",
          "Tabs are dynamically loaded from modular components.",
          "Monash Genomics and Bioinformatics Platform. ShinyCellModular: ShinyCell Package Customized by MGBP v.1 Date: Jan 2026"
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
      sc1conf_atac     = sc1conf_atac,
      sc1def_atac      = sc1def_atac,
      sc1gene_atac     = sc1gene_atac,
      sc1meta_atac     = sc1meta_atac,
      sc1fragmentpaths = sc1fragmentpaths,
      sc1annotation    = sc1annotation,
      sc1peaks         = sc1peaks,
      sc1links         = sc1links,
      markers_list = markers_list,
      assays = assays,
      dir_inputs = dir_inputs
    )

    keep <- intersect(names(args_to_pass), names(formals(srv)))
    do.call(srv, args_to_pass[keep])
  })
}


shinyApp(ui, server)
'

# I cant use sprintf to print template. error: 'fmt' length exceeds maximal format length 8192
#app_modules <- sprintf(template_app,app_title, dir_inputs, assays_str, enabled_tabs_str)

template_app <- gsub("__APP_TITLE__",    app_title,        template_app, fixed = TRUE)
template_app <- gsub("__DIR_INPUTS__",   dir_inputs,       template_app, fixed = TRUE)
template_app <- gsub("__ASSAYS__",       assays_str,       template_app, fixed = TRUE)
template_app <- gsub("__ENABLED_TABS__", enabled_tabs_str, template_app, fixed = TRUE)

app_path <- file.path(shiny.dir, "app.R")
writeLines(template_app, con = app_path)
#writeLines(app_modules, con = app_path)
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
  cat(paste0("  ", targets, " -> ", file.path(dir_prefix, targets)), sep = "\n")
  cat("\n")
  
}

invisible(app_path)
}
