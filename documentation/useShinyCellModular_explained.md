# useShinyCellModular()

`useShinyCellModular()` does not run analysis and does not compute markers, embeddings, or counts. It generates and wires a runnable ShinyCellModular application from data already prepared by `prepShinyCellModular()`.

In plain terms, it takes:
- `out_dir` — a folder containing prepared `prepShinyCellModular()` output files (`sc1conf.rds`, `sc1meta.*`, etc.)
- `shinycellmodular.dir.src` — the ShinyCellModular source directory containing `modules/`
- configuration arguments (tabs, data type, title, navbar colour, deployment)

and produces a fully functional `app.R`, with the required `modules/` copied alongside it, ready to run with `shiny::runApp()` or deploy to Posit Connect / shinyapps.io.

---

## Function signature

```r
useShinyCellModular(
    out_dir,
    shinycellmodular.dir.src = NULL,
    rsconnect.deploy         = FALSE,
    data_type                = NULL,
    enabled_tabs             = NULL,
    overwrite_modules        = FALSE,
    disable_ui_server        = TRUE,
    app_title                = NULL,
    navbar_color             = "#007BA7"
)
```

Call `useShinyCellModular(help = TRUE)` (or simply omit `out_dir`) to print the same argument summary from R.

---

## Arguments

| Argument | Default | Description |
|---|---|---|
| `out_dir` | — (required) | Directory containing prepared `prepShinyCellModular()` output files. |
| `shinycellmodular.dir.src` | `NULL` | Path to the ShinyCellModular source containing `modules/`. Auto-detected via `find.package("ShinyCellModular")` if not supplied. |
| `rsconnect.deploy` | `FALSE` | Write an rsconnect manifest for deployment (`rsconnect::writeManifest()`), then rewrite file keys to be prefixed with `out_dir`. |
| `data_type` | `NULL` | Preset tab selection — one or more of the subfolder names under `modules/` (e.g. `"RNA"`, `"RNA_ATAC"`, `"SPATIAL"`). Matched against the folders actually present on disk. |
| `enabled_tabs` | `NULL` | Character vector of tab IDs to include. Can be combined with or used instead of `data_type` (see selection rules below). |
| `overwrite_modules` | `FALSE` | Remove and replace an existing `modules/` folder in `out_dir`. |
| `disable_ui_server` | `TRUE` | Rename legacy `ui.R` / `server.R` in `out_dir` to `.bak` so Shiny doesn't prioritise them over the generated `app.R`. |
| `app_title` | `NULL` | Title shown in the app navbar. **Required** — the function stops if it's missing. |
| `navbar_color` | `"#007BA7"` | Navbar background colour (hex), also used for buttons/hover states in the generated `app.R`. |

**Return value:** invisibly returns the path to the generated `app.R`. The real output is `out_dir/app.R` plus the copied `modules/` folder (and `manifest.json` if `rsconnect.deploy = TRUE`).

---

## Step by step

### 1. Help / required arguments
If `out_dir` is missing (and `help` wasn't set), help is forced on and `helpMessage` is printed. `app_title` is mandatory — missing or `NULL` stops with an explicit error.

### 2. Resolve and validate paths
- `shinycellmodular.dir.src`: if not supplied, tries `find.package("ShinyCellModular")` and checks it has a `modules/` folder; stops with instructions if it can't be found.
- Both `out_dir` and `shinycellmodular.dir.src` are normalised with `normalizePath(mustWork = TRUE)`, so nonexistent paths fail fast with a clear error.

### 3. Determine which tabs will be included
The tab catalogue is built dynamically by scanning every file under `shinycellmodular.dir.src/modules/`, grouped by its immediate subfolder (`dirname()`), so the available `data_type` values and tab IDs always reflect what's actually on disk — nothing is hard-coded.

Selection rules:
- Neither `data_type` nor `enabled_tabs` supplied → stops with an error asking for one of them.
- Only `data_type` supplied → all tabs found under that subfolder(s) are included (`data_type` supports multiple values via `several.ok = TRUE`).
- Only `enabled_tabs` supplied → every subfolder is searched for matching tab IDs; any tab not found anywhere triggers a `warning()` listing the missing ones and is dropped.
- Both supplied → `enabled_tabs` is intersected with the tabs valid for `data_type`; any tab outside that set stops execution with the list of valid tabs for that `data_type`.

### 4. Optionally disable legacy `ui.R` / `server.R`
If `disable_ui_server = TRUE` and either file exists in `out_dir`, both are renamed to `.bak` (with a message explaining why), ensuring Shiny loads the generated `app.R`.

### 5. Copy/refresh modules and write `app.R`
Copies the module files for `enabled_tabs` into `out_dir/modules/` (replacing the folder first if `overwrite_modules = TRUE`), then builds `app.R` from an internal template string. The template:
- Sources every `.R` file in `modules/`, calling `register_tab()` per module and collecting a `tab_registry`.
- Warns if no modules registered any tabs, or if a requested tab wasn't found in the registry.
- Builds one `tabPanel` per enabled tab, appending a footer with the module's `author` / `description` / `version` / `date` / `source` / `contact` metadata.
- Wires up `server()` by passing through only the arguments each module's server function actually declares (matched against `formals()`), so modules only need to list the globals they use.
- Substitutes placeholders (`__APP_TITLE__`, `__DIR_INPUTS__`, `__ASSAYS__`, `__ENABLED_TABS__`, `__NAVBAR_COLOR__`, `__SCM_VERSION__`) with `gsub(..., fixed = TRUE)` and writes the result to `out_dir/app.R`.

### 6. Optional rsconnect manifest
If `rsconnect.deploy = TRUE`: runs `rsconnect::writeManifest(appDir = out_dir)`, then reads `manifest.json` back in, re-keys every `.rds` / `.h5` / `.parquet` entry to be prefixed with `out_dir` (so paths resolve correctly for the deployed app), and rewrites the manifest.

---

## Notes

- `data_type` and `enabled_tabs` are validated against real folders/files, so a typo in either produces an explicit list of what's actually available rather than a silent no-op.
- The generated `app.R` is self-contained — it does not call back into `useShinyCellModular()` at runtime, so `shinycellmodular.dir.src` is only needed at generation time.
- Because server arguments are matched via `formals()`, adding new global objects to the template's `args_to_pass` list is safe for existing modules — modules that don't declare the new argument simply won't receive it.
