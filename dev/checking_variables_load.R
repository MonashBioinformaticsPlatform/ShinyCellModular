scm_load_vars <- function(dir_inputs) {
  assay_dirs <- list.dirs(dir_inputs, full.names = FALSE, recursive = FALSE)
  assay_dirs <- assay_dirs[file.exists(file.path(dir_inputs, assay_dirs, "sc1conf.rds"))]

  rows <- list()
  for (assay in assay_dirs) {
    adir   <- file.path(dir_inputs, assay)
    suffix <- if (toupper(assay) == "RNA") "" else paste0("_", tolower(assay))

    entries <- list.files(adir, full.names = FALSE)
    is_dir  <- file.info(file.path(adir, entries))$isdir
    files   <- entries[!is_dir]
    subdirs <- entries[is_dir]

    for (fname in files) {
      fpath   <- file.path(adir, fname)
      base    <- tools::file_path_sans_ext(fname)
      ext     <- tolower(tools::file_ext(fname))
      varname <- paste0(base, suffix)
      rows[[length(rows) + 1]] <- data.frame(
        variable = varname, assay = assay, source_path = fpath,
        kind = "file", ext = ext, stringsAsFactors = FALSE
      )
    }

    for (sd in subdirs) {
      varname <- paste0(sd, suffix)
      rows[[length(rows) + 1]] <- data.frame(
        variable = varname, assay = assay, source_path = file.path(adir, sd),
        kind = "dir", ext = NA_character_, stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# Only inspects the ui/server functions actually passed to register_tab(),
# not every helper function defined in the file — those have their own
# unrelated parameters (chr, start, end, gene, etc.) that scm_globals never touches.
scm_module_args <- function(modules_dir) {
  files <- list.files(modules_dir, pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)
  rows  <- list()

  for (f in files) {
    exprs <- tryCatch(parse(f), error = function(e) NULL)
    if (is.null(exprs)) next

    fn_defs <- list()
    register_calls <- list()

    for (e in exprs) {
      if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
          is.call(e[[3]]) && identical(as.character(e[[3]][[1]]), "function")) {
        fn_defs[[as.character(e[[2]])]] <- e[[3]]
      }
      if (is.call(e) && identical(as.character(e[[1]]), "register_tab")) {
        register_calls[[length(register_calls) + 1]] <- e
      }
    }

    for (call in register_calls) {
      nm <- names(call)
      ui_sym     <- if ("ui" %in% nm) as.character(call[["ui"]]) else NA
      server_sym <- if ("server" %in% nm) as.character(call[["server"]]) else NA
      tab_id     <- if ("id" %in% nm) tryCatch(eval(call[["id"]]), error = function(e) NA) else NA

      for (role_name in c(ui_sym, server_sym)) {
        if (is.na(role_name) || is.null(fn_defs[[role_name]])) next
        fn <- tryCatch(eval(fn_defs[[role_name]]), error = function(e) NULL)
        if (!is.function(fn)) next
        args <- names(formals(fn))
        rows[[length(rows) + 1]] <- data.frame(
          module_file = basename(f), tab_id = tab_id, fn_name = role_name,
          arg = args, stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

# ---- Run it -----------------------------------------------------------------
dir_inputs  <- "/home/lper0012/tasks/liam.kealy/processing/github/SCmultiomics_longCOVID/data/testing/"
modules_dir <- file.path(dir_inputs, "modules")

loaded   <- scm_load_vars(dir_inputs)
mod_args <- scm_module_args(modules_dir)

ignore <- c("id", "input", "output", "session", "dir_inputs")
mod_args <- mod_args[!mod_args$arg %in% ignore, ]

strip_suffix <- function(x) sub("_(atac|multi|regulon)$", "", x)

loaded$root   <- strip_suffix(loaded$variable)
mod_args$root <- strip_suffix(mod_args$arg)

arg_by_root <- aggregate(
  arg ~ root, data = mod_args,
  FUN = function(x) paste(sort(unique(x)), collapse = " | ")
)
names(arg_by_root)[2] <- "expected_by_modules"

loaded_df <- merge(loaded, arg_by_root, by = "root", all.x = TRUE)
loaded_df$expected_by_modules[is.na(loaded_df$expected_by_modules)] <- "(not referenced by any module)"

loaded_df$match <- mapply(function(v, exp) {
  if (exp == "(not referenced by any module)") return("unused")
  opts <- strsplit(exp, " \\| ")[[1]]
  if (v %in% opts) "OK" else "MISMATCH"
}, loaded_df$variable, loaded_df$expected_by_modules)

loaded_df <- loaded_df[order(loaded_df$assay, loaded_df$variable),
                       c("variable", "assay", "kind", "expected_by_modules", "match", "source_path")]
rownames(loaded_df) <- NULL

needed_df <- data.frame(arg = sort(unique(mod_args$arg)), stringsAsFactors = FALSE)
needed_df$declared_in <- sapply(needed_df$arg, function(a) {
  paste(unique(mod_args$module_file[mod_args$arg == a]), collapse = ", ")
})
needed_df$has_loaded_match <- needed_df$arg %in% loaded_df$variable
rownames(needed_df) <- NULL

missing_df <- needed_df[!needed_df$has_loaded_match, c("arg", "declared_in")]
rownames(missing_df) <- NULL

dup_names <- loaded$variable[duplicated(loaded$variable) | duplicated(loaded$variable, fromLast = TRUE)]
collisions_df <- loaded[loaded$variable %in% unique(dup_names), c("variable", "assay", "kind", "source_path")]
collisions_df <- collisions_df[order(collisions_df$variable, collisions_df$assay), ]
rownames(collisions_df) <- NULL

cat("\n--- loaded_df: variables the loop creates, and what modules expect instead ---\n")
print(loaded_df)
cat("\n--- needed_df: every argument the registered ui/server functions declare ---\n")
print(needed_df)
cat("\n--- missing_df: registered arguments with NO matching loaded variable ---\n")
if (nrow(missing_df) == 0) cat("None.\n") else print(missing_df)
cat("\n--- collisions_df: variable names written by more than one assay folder ---\n")
if (nrow(collisions_df) == 0) cat("None.\n") else print(collisions_df)

write.table(loaded_df, "loaded_df.txt",sep="\t")
write.table(needed_df, "needed_df.txt",sep="\t")
write.table(missing_df, "missing_df.txt",sep="\t")
write.table(collisions_df, "collisions_df.txt",sep="\t")

