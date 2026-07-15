#' @importFrom yaml read_yaml write_yaml
NULL

#' Build a portable, deployable bundle from an existing ShinyCellModular app folder
#'
#' All customization and testing happens in \code{out_dir} (the folder created and
#' iterated on with \code{\link{useShinyCellModular}}). When ready to containerize, put
#' behind ShinyProxy, or serve as one of several datasets on a website engine, run
#' this function to assemble a self-contained portable bundle in \code{bundle_dir}:
#' a skeleton \code{app.R}, \code{app_config.yml}, and copies of \code{modules/}
#' and the input data. \code{out_dir} itself is never modified, it stays the dev
#' copy. \code{bundle_dir} is what gets moved to wherever it needs to live.
#'
#' \code{app_title}, \code{enabled_tabs}, \code{navbar_css}, and the generating
#' \code{ShinyCellModular} version are all read by sourcing \code{out_dir}'s
#' \code{app.R} and taking the resulting variables directly, not by parsing the
#' text. \code{assays} come from \code{out_dir}'s folder structure: every
#' top-level folder other than \code{modules/} is treated as data input. Any
#' hand-added or edited tabs in \code{app.R} and \code{modules/} are picked up
#' and carried into the bundle.
#'
#' The skeleton \code{app.R} bakes in the \code{ShinyCellModular} version that
#' generated \code{out_dir}'s \code{app.R} (\code{SKELETON_VERSION}, read from
#' that \code{app.R}, not from whatever version happens to be installed right
#' now), and compares it against \code{app_config.yml}'s \code{scm_version} at
#' startup. On mismatch, e.g. an \code{app_config.yml} paired with the wrong
#' skeleton copy, it warns but still launches, the developer decides whether
#' that mismatch actually matters.
#'
#' @param out_dir Existing, already-tested app folder created by
#'   \code{\link{useShinyCellModular}} (must contain \code{app.R}). Left untouched.
#' @param bundle_dir Destination folder for the portable bundle. Default
#'   \code{NULL}: uses \code{paste0(out_dir, "_portable")}. If it already exists
#'   and does not contain an \code{app_config.yml}, this function stops rather
#'   than risk overwriting the \code{out_dir} you are trying to convert. If it
#'   already exists and does contain an \code{app_config.yml} (an earlier
#'   portable bundle), set \code{overwrite = TRUE} to replace it.
#' @param overwrite Replace an existing \code{bundle_dir} that already contains
#'   an \code{app_config.yml}. Default \code{FALSE}.
#'
#' @return Invisibly returns \code{NULL}. Writes the bundle to \code{bundle_dir}.
#' @export
makeShinyCellModularPortable <- function(
    out_dir,
    bundle_dir = NULL,
    overwrite = FALSE
) {

  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  app_r_path <- file.path(out_dir, "app.R")

  if (!file.exists(app_r_path)) {
    stop("No app.R found in ", out_dir, ". Run useShinyCellModular() first, then build a portable bundle from it.", call. = FALSE)
  }

  if (is.null(bundle_dir)) {
    bundle_dir <- paste0(out_dir, "_portable")
  }

  if (dir.exists(bundle_dir)) {
    if (!file.exists(file.path(bundle_dir, "app_config.yml"))) {
      # this check is not gated by overwrite on purpose, an existing folder with
      # no app_config.yml is likely a ShinyCellModular dev folder, not a bundle,
      # and overwriting it would destroy the app being converted
      stop(
        "bundle_dir exists and does not contain an app_config.yml: ", bundle_dir, ". ",
        "You may be about to overwrite the ShinyCellModular app you are trying to convert into a portable bundle. ",
        "This is not recommended, please use a different bundle_dir.",
        call. = FALSE
      )
    }
    if (!isTRUE(overwrite)) {
      stop("bundle_dir already exists: ", bundle_dir, ". Set overwrite = TRUE to replace it.", call. = FALSE)
    }
    unlink(bundle_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  bundle_dir <- normalizePath(bundle_dir, mustWork = TRUE)

  existing_app <- paste(readLines(app_r_path, warn = FALSE), collapse = "\n")

  ###########################################################################
  # app_title, enabled_tabs, navbar_css, and shinyCellModularVersion are all
  # real variables in app.R, read by sourcing it rather than parsing text,
  # immune to whitespace/quoting differences that would make regex fragile.
  # app.R keeps going after these are assigned (tab building, UI, server), we
  # do not need any of that, so wrapped in tryCatch: a later error (missing
  # package, etc) is fine as long as all four variables were already assigned
  # by the time it happens
  ###########################################################################

  source_env <- new.env()
  tryCatch(
    source(app_r_path, local = source_env, echo = FALSE),
    error = function(e) {
      message("Sourcing ", app_r_path, " stopped partway (", conditionMessage(e), "). ",
              "Continuing, this is fine as long as app_title/navbar_css/shinyCellModularVersion/enabled_tabs were already assigned by then.")
    }
  )

  if (is.null(source_env$app_title)) {
    stop("Could not read app_title by sourcing ", app_r_path, ".", call. = FALSE)
  }
  app_title <- source_env$app_title

  if (is.null(source_env$navbar_css)) {
    stop("Could not read navbar_css by sourcing ", app_r_path, ".", call. = FALSE)
  }
  navbar_css <- source_env$navbar_css

  if (is.null(source_env$shinyCellModularVersion)) {
    scm_version <- as.character(utils::packageVersion("ShinyCellModular"))
    warning("Could not read shinyCellModularVersion by sourcing ", app_r_path,
            ". Falling back to the currently installed package version (",
            scm_version, "), which may not be the version that actually generated this app.R.", call. = FALSE)
  } else {
    scm_version <- as.character(source_env$shinyCellModularVersion)
  }

  if (is.null(source_env$enabled_tabs)) {
    stop("Could not read enabled_tabs by sourcing ", app_r_path, ".", call. = FALSE)
  }
  enabled_tabs <- source_env$enabled_tabs

  ###########################################################################
  # assays are not read from app.R at all, out_dir's own folder structure is
  # the source of truth: every top-level folder other than modules/ is data
  # input, whatever it's named, same "everything not modules/ is data" rule
  ###########################################################################

  top_level_dirs <- list.dirs(out_dir, recursive = FALSE, full.names = FALSE)
  assays <- setdiff(top_level_dirs, "modules")
  if (length(assays) == 0) {
    stop("Could not find any data folders (folders other than 'modules') inside out_dir: ", out_dir, call. = FALSE)
  }

  ###########################################################################
  # copy modules/ verbatim, this is where any hand-added/edited tab files live
  ###########################################################################

  src_modules_dir <- file.path(out_dir, "modules")
  if (!dir.exists(src_modules_dir)) {
    stop("Could not find 'modules' folder in out_dir: ", src_modules_dir, call. = FALSE)
  }
  dst_modules_dir <- file.path(bundle_dir, "modules")
  dir.create(dst_modules_dir, recursive = TRUE, showWarnings = FALSE)
  module_files <- list.files(src_modules_dir, full.names = TRUE)
  ok <- file.copy(module_files, dst_modules_dir, recursive = FALSE)
  if (!all(ok)) {
    stop("Failed to copy the following module files into ", dst_modules_dir, ":\n  ",
         paste(basename(module_files[!ok]), collapse = "\n  "))
  }
  message("Copied ", sum(ok), " module(s) into: ", dst_modules_dir)

  ###########################################################################
  # copy the input data, "the bundle of inputs"
  ###########################################################################

  for (assay_dir in assays) {
    src <- file.path(out_dir, assay_dir)
    dst <- file.path(bundle_dir, assay_dir)
    dir.create(dst, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(list.files(src, full.names = TRUE), dst, recursive = TRUE)
    if (!all(ok)) {
      stop("Failed to copy some files from ", src, " into ", dst, call. = FALSE)
    }
    message("Copied ", assay_dir, " data into: ", dst)
  }

  ###########################################################################
  # write app_config.yml
  ###########################################################################

  app_config <- list(
    app_title    = app_title,
    dir_inputs   = ".",
    assays       = as.list(assays),
    enabled_tabs = as.list(enabled_tabs),
    navbar_css   = navbar_css,
    scm_version  = scm_version
  )

  config_path <- file.path(bundle_dir, "app_config.yml")
  yaml::write_yaml(app_config, config_path)
  message("Wrote app_config.yml to: ", config_path)

  ###########################################################################
  # turn the ACTUAL existing app.R text into the skeleton, in place, textually.
  # this part is unavoidably text-editing, we are inserting new code and
  # rewriting specific lines, sourcing above only applies to reading values,
  # not to modifying the file itself
  ###########################################################################

  skeleton <- existing_app

  config_load_block <- paste0(
    '# app.R does not change between datasets in portable mode\n',
    '# every dataset just needs its own app_config.yml (and data) sitting next to a copy of this file\n',
    '# SCMODULAR_DATA_DIR still works exactly as before, it now also tells app.R where the bundle (data + app_config.yml) lives\n',
    'get_app_dir <- function() {\n',
    '  ofile <- tryCatch({sys.frame(1)$ofile}, error = function(e) {return(NULL)})\n',
    '  if (!is.null(ofile) && is.character(ofile) && nzchar(ofile)) {\n',
    '    return(normalizePath(dirname(ofile), winslash = "/", mustWork = TRUE))\n',
    '  }\n',
    '\n',
    '  args <- commandArgs(trailingOnly = FALSE)\n',
    '  file_arg <- grep("^--file=", args, value = TRUE)\n',
    '  if (length(file_arg) == 1) {\n',
    '    f <- sub("^--file=", "", file_arg)\n',
    '    if (nzchar(f) && file.exists(f)) {\n',
    '      return(normalizePath(dirname(f), winslash = "/", mustWork = TRUE))\n',
    '    }\n',
    '  }\n',
    '\n',
    '  normalizePath(getwd(), winslash = "/", mustWork = TRUE)\n',
    '}\n',
    '\n',
    'bundle_dir <- Sys.getenv("SCMODULAR_DATA_DIR", unset = get_app_dir())\n',
    'config_path <- file.path(bundle_dir, "app_config.yml")\n',
    '\n',
    'if (!file.exists(config_path)) {\n',
    '  stop("app_config.yml not found: ", config_path, ". This skeleton app.R needs an app_config.yml next to it (or point SCMODULAR_DATA_DIR at a bundle that has one).", call. = FALSE)\n',
    '}\n',
    '\n',
    'app_config <- yaml::read_yaml(config_path)\n',
    '\n',
    '# baked in at bundle-build time by makeShinyCellModularPortable(), warns (does not\n',
    '# stop) if this bundle ends up paired with a mismatched skeleton copy\n',
    'SKELETON_VERSION <- "', scm_version, '"\n',
    'if (!is.null(app_config$scm_version) && !identical(as.character(app_config$scm_version), SKELETON_VERSION)) {\n',
    '  warning(\n',
    '    "Skeleton app.R version (", SKELETON_VERSION, ") does not match app_config.yml scm_version (",\n',
    '    app_config$scm_version, "). This bundle may have been built with a different ShinyCellModular version.",\n',
    '    call. = FALSE\n',
    '  )\n',
    '}\n',
    '\n',
    'app_title    <- app_config$app_title\n',
    'navbar_css   <- app_config$navbar_css\n',
    'shinyCellModularVersion <- app_config$scm_version\n',
    '\n',
    'dir_inputs <- file.path(bundle_dir, app_config$dir_inputs)'
  )

  # app_title / navbar_css / shinyCellModularVersion / dir_inputs assignment
  # block -> read from app_config instead
  skeleton_new <- sub(
    'app_title <- "[^"]*"\n\nnavbar_css <- "[^"]*"\n\nshinyCellModularVersion <- "[^"]*"\n\ndir_inputs <- "[^"]*"',
    config_load_block,
    skeleton
  )
  if (identical(skeleton_new, skeleton)) {
    stop("Could not find the app_title/navbar_css/shinyCellModularVersion/dir_inputs block in ", app_r_path, " to make portable. ",
         "Was this app.R already made portable, or hand-edited beyond recognition?", call. = FALSE)
  }
  skeleton <- skeleton_new

  # assays <- c(...) -> read from app_config, keep whatever trailing comment was there
  skeleton <- sub(
    '(assays <- )c\\([^)]*\\)([^\n]*)',
    "\\1unlist(app_config$assays)\\2",
    skeleton
  )

  # enabled_tabs <- c(...) -> read from app_config, keep whatever trailing comment was there
  skeleton <- sub(
    '(enabled_tabs <- )c\\([^)]*\\)([^\n]*)',
    "\\1unlist(app_config$enabled_tabs)\\2",
    skeleton
  )

  # the footer already reads `shinyCellModularVersion` as a variable (see app.R's
  # template), and that variable is reassigned from app_config above, so no
  # separate footer text substitution is needed here

  app_r_out <- file.path(bundle_dir, "app.R")
  writeLines(skeleton, app_r_out)
  message("Wrote portable app.R to: ", app_r_out)

  message("Portable bundle ready at: ", bundle_dir, ". out_dir was not modified.")

  invisible(NULL)
}
