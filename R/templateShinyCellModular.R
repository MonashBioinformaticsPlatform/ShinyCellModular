#' Create a template file for a new ShinyCellModular module/tab
#'
#' Writes a skeleton .R file following the standard module structure
#' (Functions / UI / Server / Registration) used across ShinyCellModular
#' modules, so a new tab can be built by filling in the placeholders instead
#' of copy-pasting an existing module by hand.
#'
#' @param id Tab id — used as the filename (sans extension) and as the prefix
#'   for the generated function names (\code{<id>_ui}, \code{<id>_server}).
#'   Must be unique across \code{modules/}. Required.
#' @param title Human readable title shown on the tab. Required.
#' @param data_type Which \code{modules/} subfolder to write into: \code{'RNA'},
#'   \code{'RNA_ATAC'} or \code{'SPATIAL'}. Default: \code{'RNA'}.
#' @param multi Write a \code{_multi} variant. \code{_multi} is appended to
#'   \code{id} for the filename, function names and registered tab id, per the
#'   existing multi-dataset naming convention. Default: \code{FALSE}.
#' @param description Short description used in \code{register_tab()}. Default:
#'   \code{"TODO: describe this tab"}.
#' @param author Author name used in \code{register_tab()}. Required.
#' @param contact Contact email used in \code{register_tab()}. Required.
#' @param shinycellmodular.dir.src Path to the ShinyCellModular source directory
#'   containing \code{modules/}. Default: \code{system.file('', package = 'ShinyCellModular')}.
#' @param overwrite Overwrite an existing file with the same name. Default: \code{FALSE}.
#'
#' @return Invisibly returns the path to the new file.
#'
#' @examples
#' \dontrun{
#' templateShinyCellModular(
#'   id      = "my_new_tab",
#'   title   = "My New Tab",
#'   data_type = "RNA",
#'   author  = "Your Name",
#'   contact = "your.email@monash.edu"
#' )
#' }
#'
#' @export
templateShinyCellModular <- function(
    id, # tab id, becomes the filename and the <id>_ui / <id>_server prefix
    title, # human readable title shown on the tab
    data_type = "RNA", # which modules/ subfolder to write into
    multi = FALSE, # write a _multi variant
    description = "TODO: describe this tab", # short description for register_tab()
    author = NULL, # author for register_tab() — required
    contact = NULL, # contact for register_tab() — required
    shinycellmodular.dir.src = NULL, # path to ShinyCellModular source containing modules/
    overwrite = FALSE # overwrite existing file with the same name

) {

  ###########################################################################
  # Help
  ###########################################################################

  if (missing(id) || is.null(id) || !nzchar(id)) {
    stop("id is missing. Please provide a tab id, e.g. templateShinyCellModular(id = 'my_new_tab', title = 'My New Tab').")
  }
  if (missing(title) || is.null(title) || !nzchar(title)) {
    stop("title is missing. Please provide a human readable title, e.g. title = 'My New Tab'.")
  }
  if (missing(author) || is.null(author) || !nzchar(author)) {
    stop("author is missing. Please provide your name, e.g. author = 'Your Name'.")
  }
  if (missing(contact) || is.null(contact) || !nzchar(contact)) {
    stop("contact is missing. Please provide your contact email, e.g. contact = 'your.email@monash.edu'.")
  }

  data_type <- match.arg(data_type, choices = c("RNA", "RNA_ATAC", "SPATIAL"))

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
  shinycellmodular.dir.src <- normalizePath(shinycellmodular.dir.src, mustWork = TRUE)

  # _multi is baked into the id itself, matching the existing naming convention
  # (e.g. bubble_heatmap.R -> bubble_heatmap_multi.R, id = "bubble_heatmap_multi")
  module_id <- if (isTRUE(multi)) paste0(id, "_multi") else id

  dst_dir  <- file.path(shinycellmodular.dir.src, "modules", data_type)
  dst_file <- file.path(dst_dir, paste0(module_id, ".R"))

  if (!dir.exists(dst_dir)) {
    stop("Target modules folder not found: ", dst_dir)
  }

  if (file.exists(dst_file) && !isTRUE(overwrite)) {
    stop(
      "A module file already exists at: ", dst_file, "\n",
      "Set overwrite = TRUE if you want to replace it.",
      call. = FALSE
    )
  }

  ###########################################################################
  # Template
  ###########################################################################

  template_module <- '
# __DESCRIPTION__
# id     = "__ID__"
# title  = "__TITLE__"

############################################### Functions ############################################

# helper / plot functions go here
# function names should be unique enough not to collide across modules
# prefer sc_* prefix for shared-looking helpers defined locally

############################################### UI ###################################################

__ID___ui <- function(id, sc1conf, sc1def) {
  ns <- NS(id)
  # shiny UI — use ns() for all input/output IDs
  # sidebarLayout + sidebarPanel / mainPanel is the standard pattern
  tabPanel(
    "__TITLE__"
  )
}

############################################### Server ##############################################

__ID___server <- function(id, sc1conf, sc1def, sc1meta, sc1gene, inpH5, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    # reactive logic, renderPlot, renderDT, downloadHandler …
    # for _multi tabs, sc1conf / sc1meta / sc1gene / inpH5 are named lists keyed by dataset index
  })
}

############################################### Registration #########################################

register_tab(
  id          = "__ID__",
  title       = "__TITLE__",
  ui          = __ID___ui,
  server      = __ID___server,
  author      = "__AUTHOR__",
  description = "__DESCRIPTION__",
  version     = "1.0",
  date        = "__DATE__",
  source      = "MGBP custom",
  contact     = "__CONTACT__"
)
'

  template_module <- gsub("__ID__",          module_id,   template_module, fixed = TRUE)
  template_module <- gsub("__TITLE__",       title,       template_module, fixed = TRUE)
  template_module <- gsub("__AUTHOR__",      author,      template_module, fixed = TRUE)
  template_module <- gsub("__CONTACT__",     contact,     template_module, fixed = TRUE)
  template_module <- gsub("__DESCRIPTION__", description, template_module, fixed = TRUE)
  template_module <- gsub("__DATE__",        format(Sys.Date(), "%b %Y"), template_module, fixed = TRUE)

  writeLines(template_module, con = dst_file)
  message("Wrote module template to: ", dst_file)

  invisible(dst_file)
}
