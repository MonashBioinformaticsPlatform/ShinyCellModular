#' @importFrom jsonlite toJSON
#' @importFrom stats na.omit
#' @importFrom utils install.packages
NULL
# packages used inside generated app.R template — referenced here to satisfy R CMD check
if (FALSE) {
  jsonlite::toJSON
  rsconnect::deployApp
  library(jsonlite)
  library(rsconnect)
}
# not exported, only used by roxygen at document() time to build the
# data_type param doc from whatever folders actually exist under modules/
#' @noRd
scm_data_type_doc <- function() {
  root <- system.file("modules", package = "ShinyCellModular")
  if (root == "" || !dir.exists(root)) {
    return("@param data_type Preset tab selection. One or more subfolder names under `modules/`.")
  }
  types <- list.dirs(root, full.names = FALSE, recursive = FALSE)
  paste0(
    "@param data_type Preset tab selection: one or more of ",
    paste0("`", types, "`", collapse = ", "),
    ". Matched against the folders actually present on disk under `modules/`."
  )
}

#' Generate a modular ShinyCellModular Shiny app
#'
#' Takes the output files from \code{\link{prepShinyCellModular}} and writes a
#' ready-to-run \code{app.R} with the selected module tabs into \code{out_dir}.
#'
#' @param out_dir Directory containing prepared prepShinyCellModular output files
#'   (\code{sc1conf.rds}, \code{sc1meta.rds}, etc.).
#' @param shinycellmodular.dir.src Path to the ShinyCellModular source directory
#'   containing \code{modules/}. Default: \code{system.file('', package = 'ShinyCellModular')}.
#' @param rsconnect.deploy Write an rsconnect manifest for deployment. Default: \code{FALSE}.
#' @eval scm_data_type_doc()
#' @param enabled_tabs Character vector of tab IDs to include. Overrides \code{data_type}
#'   presets. Default: \code{NULL} (uses all tabs for \code{data_type}).
#' @param overwrite_modules Remove and replace existing \code{modules/} folder. Default: \code{FALSE}.
#' @param disable_ui_server Rename legacy \code{ui.R} and \code{server.R} to \code{.bak}. Default: \code{TRUE}.
#' @param app_title Title shown in the app navbar. Required.
#' @param navbar_color Navbar background colour (hex), also used for buttons/hover
#'   states in the generated \code{app.R}. Default: \code{"#007BA7"}.
#'
#' @return Invisibly returns \code{NULL}. Writes \code{app.R} and \code{modules/} to \code{out_dir}.
#'
#' @section Step by step:
#' \enumerate{
#'   \item \strong{Help / required arguments.} If \code{out_dir} is missing (and
#'     \code{help} wasn't set), help is forced on and the help message is printed.
#'     \code{app_title} is mandatory, missing or \code{NULL} stops with an explicit
#'     error.
#'   \item \strong{Resolve and validate paths.} If \code{shinycellmodular.dir.src}
#'     isn't supplied, tries \code{find.package("ShinyCellModular")} and checks it has
#'     a \code{modules/} folder, stopping with instructions if it can't be found. Both
#'     \code{out_dir} and \code{shinycellmodular.dir.src} are normalised with
#'     \code{normalizePath(mustWork = TRUE)}, so nonexistent paths fail fast with a
#'     clear error.
#'   \item \strong{Determine which tabs will be included.} The tab catalogue is built
#'     dynamically by scanning every file under
#'     \code{shinycellmodular.dir.src/modules/}, grouped by its immediate subfolder, so
#'     the available \code{data_type} values and tab IDs always reflect what's
#'     actually on disk, nothing is hard-coded. Selection rules:
#'     \itemize{
#'       \item Neither \code{data_type} nor \code{enabled_tabs} supplied: stops with an
#'         error asking for one of them.
#'       \item Only \code{data_type} supplied: all tabs found under that subfolder(s)
#'         are included (\code{data_type} supports multiple values via
#'         \code{several.ok = TRUE}).
#'       \item Only \code{enabled_tabs} supplied: every subfolder is searched for
#'         matching tab IDs, any tab not found anywhere triggers a \code{warning()}
#'         listing the missing ones and is dropped.
#'       \item Both supplied: \code{enabled_tabs} is intersected with the tabs valid
#'         for \code{data_type}, any tab outside that set stops execution with the list
#'         of valid tabs for that \code{data_type}.
#'     }
#'   \item \strong{Optionally disable legacy \code{ui.R} / \code{server.R}.} If
#'     \code{disable_ui_server = TRUE} and either file exists in \code{out_dir}, both
#'     are renamed to \code{.bak} (with a message explaining why), ensuring Shiny loads
#'     the generated \code{app.R}.
#'   \item \strong{Copy/refresh modules and write \code{app.R}.} Copies the module
#'     files for \code{enabled_tabs} into \code{out_dir/modules/} (replacing the folder
#'     first if \code{overwrite_modules = TRUE}), then builds \code{app.R} from an
#'     internal template string. The template:
#'     \itemize{
#'       \item Sources every \code{.R} file in \code{modules/}, calling
#'         \code{register_tab()} per module and collecting a \code{tab_registry}.
#'       \item Warns if no modules registered any tabs, or if a requested tab wasn't
#'         found in the registry.
#'       \item Builds one \code{tabPanel} per enabled tab, appending a footer with the
#'         module's \code{author} / \code{description} / \code{version} / \code{date} /
#'         \code{source} / \code{contact} metadata.
#'       \item Wires up \code{server()} by passing through only the arguments each
#'         module's server function actually declares (matched against
#'         \code{formals()}), so modules only need to list the globals they use.
#'       \item Substitutes placeholders (\code{__APP_TITLE__}, \code{__DIR_INPUTS__},
#'         \code{__ASSAYS__}, \code{__ENABLED_TABS__}, \code{__NAVBAR_COLOR__},
#'         \code{__SCM_VERSION__}) with \code{gsub(..., fixed = TRUE)} and writes the
#'         result to \code{out_dir/app.R}.
#'     }
#'   \item \strong{Optional rsconnect manifest.} If \code{rsconnect.deploy = TRUE}:
#'     runs \code{rsconnect::writeManifest(appDir = out_dir)}, then reads
#'     \code{manifest.json} back in, re-keys every \code{.rds} / \code{.h5} /
#'     \code{.parquet} entry to be prefixed with \code{out_dir} (so paths resolve
#'     correctly for the deployed app), and rewrites the manifest.
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item \code{data_type} and \code{enabled_tabs} are validated against real
#'     folders/files, so a typo in either produces an explicit list of what's actually
#'     available rather than a silent no-op.
#'   \item The generated \code{app.R} is self-contained, it does not call back into
#'     \code{useShinyCellModular()} at runtime, so \code{shinycellmodular.dir.src} is
#'     only needed at generation time.
#'   \item Because server arguments are matched via \code{formals()}, adding new global
#'     objects to the template's \code{args_to_pass} list is safe for existing modules,
#'     modules that don't declare the new argument simply won't receive it.
#' }
#'
#' @examples
#' \dontrun{
#' useShinyCellModular(
#'   out_dir  = "my_app_files/",
#'   data_type  = "RNA",
#'   app_title  = "My scRNA-seq App"
#' )
#' }
#'
#' @export
useShinyCellModular <- function(
    out_dir, # files from shinycell are
    shinycellmodular.dir.src = NULL, # modules where shinycellmodular 
    rsconnect.deploy = FALSE, # do you want to publish in rsconnect
    data_type = NULL, # what predetermine tabs you want
    enabled_tabs = NULL, # what tabs you want
    overwrite_modules = FALSE, # overwrite modules
    disable_ui_server = TRUE, # this disables the existing ui.R and server.r
    app_title=NULL,
    navbar_color = "#00BDC9" # navbar background colour
    
) {
  
  ###########################################################################
  # Help
  ###########################################################################
  
  if (missing(out_dir) && !isTRUE(help)) {
    help <- TRUE
  }
  if (isTRUE(help)) {
    tools::Rd2txt(tools::Rd_db("ShinyCellModular")[["useShinyCellModular.Rd"]])
    return(invisible(NULL))
  }
  
  if (is.null(shinycellmodular.dir.src) || !nzchar(shinycellmodular.dir.src)) {
    # try installed package first, then fall back to find.package() which also works for load_all()
    pkg_dir <- tryCatch(find.package("ShinyCellModular"), error = function(e) "")
    if (nzchar(pkg_dir) && dir.exists(file.path(pkg_dir, "modules"))) {
      shinycellmodular.dir.src <- pkg_dir
    } else {
      stop(
        "Could not auto-detect shinycellmodular.dir.src. ",
        "Please pass it explicitly, e.g.:\n",
        "  shinycellmodular.dir.src = system.file('', package = 'ShinyCellModular')\n",
        "or the path to your local modules of ShinyCellModular.",
        call. = FALSE
      )
    }
    message("Auto-detected shinycellmodular.dir.src: ", shinycellmodular.dir.src)
  }
  # missing() above does not catch out_dir = "", check separately
  if (!nzchar(out_dir)) {
    stop("out_dir is empty. Please provide the directory containing your prepared prepShinyCellModular files.")
  }
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  shinycellmodular.dir.src <- normalizePath(shinycellmodular.dir.src, mustWork = TRUE)
  
  message("ShinyCellModular app generation starting")
  message("Target app directory: ", out_dir)
  message("Source ShinyCellModular directory: ", shinycellmodular.dir.src)
  

  
  if (missing(app_title) || is.null(app_title)) {
    stop("App title missing. We are not launching anonymous software today. Please provide a title using app_title='...'.")
  }

#### Handling the different data types and enable tabs 
  
  # treat NULL, "", and character(0) as empty
  is_empty <- function(x) {
    is.null(x) || length(x) == 0 || (is.character(x) && all(trimws(x) == ""))
  }
  data_type_provided <- !is_empty(data_type)
  tabs_provided      <- !is_empty(enabled_tabs)
  
  if (!data_type_provided && !tabs_provided) {
    stop("You must provide either data_type or enabled_tabs. For example data_type can be 'RNA' or RNA_ATAC or SPATIAL. What type of assay do you have?")
  }
  
  # tab catalogue: every tab ID per data_type — single source of truth
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
  
  if (data_type_provided) {
    data_type <- match.arg(data_type, choices = names(all_tabs_by_type),several.ok = TRUE)
  }
  
  if (data_type_provided && !tabs_provided) {
    # all tabs for that data_type
    enabled_tabs <- unique(unlist(lapply(data_type, function(dt) all_tabs_by_type[[dt]]$tab_id)))
    
  } else if (!data_type_provided && tabs_provided) {
    # search all folders for requested tabs
    found <- unlist(lapply(names(all_tabs_by_type), function(dt) {
      intersect(enabled_tabs, all_tabs_by_type[[dt]]$tab_id)
    }))
    missing_tabs <- setdiff(enabled_tabs, found)
    if (length(missing_tabs) > 0)
      warning("The following tabs were not found in any data_type folder and will be skipped: ",
              paste(missing_tabs, collapse = ", "), call. = FALSE)
    enabled_tabs <- found
    
  } else {
    # both provided — intersection of data_type tabs and requested tabs
    allowed_tabs <- all_tabs_by_type[[data_type]]$tab_id
    bad_tabs     <- setdiff(enabled_tabs, allowed_tabs)
    if (length(bad_tabs) > 0) {
      stop(
        "The following tabs are not valid for data_type = '", data_type, "':\n",
        "  ", paste(bad_tabs, collapse = ", "), "\n\n",
        "  allowed tabs for '", data_type, "':\n",
        "  ", paste(allowed_tabs, collapse = ", "), "\n\n",
        "To use tabs from a different data_type, change data_type accordingly.",
        call. = FALSE
      )
    }
    enabled_tabs <- intersect(allowed_tabs, enabled_tabs)
  }
  
  message("Enabled tabs : ", paste(enabled_tabs, collapse = ", "))
#####
  
  
  if (isTRUE(disable_ui_server)) {
    ui_r <- file.path(out_dir, "ui.R")
    server_r <- file.path(out_dir, "server.R")
    
    if (file.exists(ui_r) || file.exists(server_r)) {
      message(
        "ui.R and or server.R detected in the app directory. ",
        "Shiny will prioritise these files over app.R. ",
        "To ensure the modular ShinyCellModular app is used, ",
        "ui.R and server.R will be disabled by renaming them. ",
        "Backup files with extension .bak will be created."
      )
    }
    
    if (file.exists(ui_r)) {
      file.rename(ui_r, file.path(out_dir, "ui.R.bak"))
      message("Renamed ui.R to ui.R.bak")
    }
    
    if (file.exists(server_r)) {
      file.rename(server_r, file.path(out_dir, "server.R.bak"))
      message("Renamed server.R to server.R.bak")
    }
  }
  
  
  #idx_tabs<-which(all_tabs_by_type[[data_type]]$tab_id %in% enabled_tabs)
  
  # build src_modules — tabs may come from multiple folders when only enabled_tabs is provided
  src_modules <- unlist(lapply(names(all_tabs_by_type), function(dt) {
    idx <- which(all_tabs_by_type[[dt]]$tab_id %in% enabled_tabs)
    if (length(idx) == 0) return(character(0))
    file.path(shinycellmodular.dir.src, "modules", dt, all_tabs_by_type[[dt]]$filename[idx])
  }))
  
  # copy exact requested tabs only, avoid copying excessive amount of files if this app grow
  #src_modules <- file.path(shinycellmodular.dir.src, "modules/",data_type,all_tabs_by_type[[data_type]]$filename[idx_tabs])
  src_modules_dir <- file.path(shinycellmodular.dir.src, "modules/")
  dst_modules <- file.path(out_dir, "modules/")
  
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
    message("You did not ask me to overwrite the modules folder. If you want me to rewrite them set overwrite_modules = TRUE. For now using existing modules/ folder in: ", dst_modules)
  }
  
  #dir_inputs <- out_dir
  dir_inputs <- Sys.getenv("SCMODULAR_DATA_DIR", unset = out_dir)
  
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
    axis.ticks =  element_line(colour = "black", linewidth = base_size / 20),
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

if (file.exists(file.path(dir_inputs,"RNA"))) {
        rna_dir         <- file.path(dir_inputs, "RNA")
        sc1conf <- tryCatch({readRDS(file.path(rna_dir, "sc1conf.rds"))},error = function(e) {return(NULL)})
        sc1def  <- tryCatch({readRDS(file.path(rna_dir, "sc1def.rds"))},error = function(e) {return(NULL)})
        sc1gene <- tryCatch({readRDS(file.path(rna_dir, "sc1gene.rds"))},error = function(e) {return(NULL)})
        sc1meta <- tryCatch({readRDS(file.path(rna_dir, "sc1meta.rds"))},error = function(e) {return(NULL)})
        markers_list <- tryCatch({file.path(rna_dir, "markergenes_lists.parquet")},error = function(e) {return(NULL)})

}else { sc1conf     <- NULL
       sc1def      <- NULL
       sc1gene     <- NULL
       sc1meta     <- NULL
       markers_list <- NULL
        }

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
  sc1atactiles     <- tryCatch({ d <- file.path(atac_dir, "tiles")
                         if (dir.exists(d)) d else NULL
                       }, error = function(e) {return(NULL)})
}else { sc1conf_atac     <- NULL
          sc1def_atac      <- NULL
          sc1gene_atac     <- NULL
          sc1meta_atac     <- NULL
          sc1fragmentpaths <- NULL
          sc1annotation    <- NULL
          sc1peaks         <- NULL
          sc1links         <- NULL
          sc1atactiles     <- NULL
        }


assays <- __ASSAYS__ # still unclear if I am using this for anything
assays_vec <- unique(sc1conf$assay)

tab_registry <- list()

register_tab <- function(id, title, ui, server,
                         author = NULL, description = NULL, version = NULL,
                         date = NULL, source = NULL, contact = NULL) {
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
    title       = title,
    ui          = ui,
    server      = server,
    # --- tab metadata (shown in per-tab footer) ---
    author      = author,
    description = description,
    version     = version,
    date        = date,
    source      = source,
    contact     = contact
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

  meta <- tab_registry[[k]]

  # build per-tab footer from registry metadata — always shown, fields omitted if NULL
  .fmt_field <- function(label, val) {
    if (!is.null(val) && nzchar(val)) tags$span(style = "margin-right: 18px;", tags$b(paste0(label, ": ")), val)
  }
  tab_footer <- tags$div(
    style = paste(
      "margin-top: 32px; padding: 8px 14px; border-top: 1px solid #ddd;",
      "background: #f8f8f8; font-size: 82%; color: #888; line-height: 1.8;"
    ),
    tags$b("Tab info:"), tags$br(),
    .fmt_field("Author",      meta$author),
    .fmt_field("Description", meta$description),
    .fmt_field("Version",     meta$version),
    .fmt_field("Date",        meta$date),
    .fmt_field("Source",      meta$source),
    .fmt_field("Contact",     meta$contact)
  )

  tabPanel(
    meta$title,
    meta$ui(id = k, sc1conf = sc1conf, sc1def = sc1def),
    tab_footer
  )
})



 ui <- fluidPage( 
      tags$head(
        tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}")),
        tags$style(HTML(".navbar-default{background-color:__NAVBAR_COLOR__;border-color:__NAVBAR_COLOR__;} .navbar{min-height:36px;font-family:Helvetica,Arial,sans-serif;} .navbar-default .navbar-nav>li>a{color:#fff;padding-top:8px;padding-bottom:8px;font-size:13px;font-weight:bold;letter-spacing:0.3px;border-radius:4px;transition:background-color 0.15s ease;} .navbar-default .navbar-nav>li>a:hover{background-color:rgba(255,255,255,0.15);} .navbar-default .navbar-brand{color:#fff;padding-top:8px;padding-bottom:8px;font-size:13px;font-weight:bold;letter-spacing:0.3px;border-right:1px solid rgba(255,255,255,0.3);margin-right:4px;} .navbar-collapse{padding-top:0;padding-bottom:0;} .btn-default{background-color:__NAVBAR_COLOR__;border-color:__NAVBAR_COLOR__;color:#fff;border-radius:6px;font-weight:bold;transition:filter 0.15s ease;} .btn-default:hover,.btn-default:focus,.btn-default:active{background-color:__NAVBAR_COLOR__;border-color:__NAVBAR_COLOR__;color:#fff;filter:brightness(0.92);}"))

    ),
      do.call(navbarPage, c(list(title = app_title), tab_panels)),
      tags$hr(),
tags$p(
  style = "font-size: 90%; color: #666;",
  em(
    "This application was generated using ",
    tags$a(
      "ShinyCellModular",
      href   = "https://github.com/MonashBioinformaticsPlatform/ShinyCellModular",
      target = "_blank"
    ),
    paste0(
      " v__SCM_VERSION__",
      ": a modular Shiny framework for single-cell data exploration developed by the ",
      "Monash Genomics and Bioinformatics Platform (MGBP), extending the ShinyCell package"
    )
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
      sc1atactiles     = sc1atactiles,
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
shinyCellModularVersion<-packageVersion("ShinyCellModular")

template_app <- gsub("__APP_TITLE__",    app_title,        template_app, fixed = TRUE)
template_app <- gsub("__DIR_INPUTS__",   dir_inputs,       template_app, fixed = TRUE)
template_app <- gsub("__ASSAYS__",       assays_str,       template_app, fixed = TRUE)
template_app <- gsub("__ENABLED_TABS__", enabled_tabs_str, template_app, fixed = TRUE)
template_app <- gsub("__NAVBAR_COLOR__", navbar_color, template_app, fixed = TRUE)
template_app <- gsub("__SCM_VERSION__", shinyCellModularVersion, template_app, fixed = TRUE)



app_path <- file.path(out_dir, "app.R")
writeLines(template_app, con = app_path)
#writeLines(app_modules, con = app_path)
message("Wrote app.R to: ", app_path)


if (isTRUE(rsconnect.deploy)) {
  library(rsconnect)
  library(jsonlite)
  
  rsconnect::writeManifest(appDir = out_dir)
  message("Wrote rsconnect manifest in: ", out_dir)
  
  dir_prefix <- out_dir
  manifest_path <- file.path(out_dir, "manifest.json")
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
