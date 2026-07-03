---
name: shinycellmodular
description: "Development skill for the ShinyCellModular R package, a modular Shiny app framework for single-cell (scRNA-seq, scATAC-seq, multiome, spatial) data built on top of ShinyCell. Use this skill whenever the user is adding a new module/tab, modifying an existing tab, adding a feature to a tab, modifying prepShinyCellModular, modifying useShinyCellModular, debugging a module, or adding a new data_type. Also use when the user asks how the project is structured, how register_tab works, how modules are loaded, or how prep/use functions fit together, even if they don't say 'ShinyCellModular' explicitly."
---

# ShinyCellModular Development Skill

## Project Overview

ShinyCellModular wraps [ShinyCell](https://github.com/SGDDNB/ShinyCell) output into a
fully modular Shiny app. Each visualisation tab is a self-contained R file that registers
itself at source-time via `register_tab()`. Two top-level functions drive the workflow:

| Function | File | Role |
|---|---|---|
| `prepShinyCellModular()` | `prepShinyCellModular.R` | Pre-processes a Seurat object → files on disk |
| `useShinyCellModular()` | `useShinyCellModular.R` | Writes `app.R` + copies `modules/` to the target Shiny dir |

---

## Repository Layout

```
ShinyCellModular/
├── R/
│   ├── prepShinyCellModular.R   # Seurat → ShinyCell files + extras (markers, ATAC, 3D UMAP …)
│   └── useShinyCellModular.R    # app.R generator; copies modules/; writes manifest.json
├── inst/
│   └── modules/
│       ├── RNA/                 # single-assay RNA tabs
│       ├── MULTI/            # multiome tabs (RNA + ATAC in same object)
│       └── SPATIAL/             # spatial tabs
└── DESCRIPTION / NAMESPACE / …
```

Each file in `modules/<data_type>/` is one tab. The filename (sans extension) is the tab's
`id`. `useShinyCellModular()` discovers available tabs by listing that directory at run time —
no hard-coded catalogue.

---

## Adding a Feature to an Existing Tab (most common request)

This is the task this skill is most often used for. Follow this order:

1. **Find the one file.** Open `inst/modules/<data_type>/<tab_id>.R`. That file is the
   entire tab — UI, server, and registration. You should not need to open, or change,
   any other module file to add a feature to this one.
2. **Read the whole file before editing.** Note the existing helper functions, the UI
   layout pattern, and which server arguments (see table below) the tab already uses.
3. **Make the minimum-diff change.** Add the new input/output, the new reactive logic,
   and any new helper function — in the same style as what's already there. Do not
   reformat, rename, or refactor surrounding code that wasn't asked about.
4. **If the feature needs extra files** beyond what prep already writes, they just need to
   be dropped somewhere the tab can find them (alongside `sc1conf.rds`, `sc1counts.h5`,
   etc. in the prep output directory), and read by the module's server function.
   Automating that generation inside `prepShinyCellModular.R` is something to consider
   later, not a requirement for adding the feature.
5. **Never touch `useShinyCellModular.R`** for a tab-level feature. It doesn't need to
   know about individual tabs' internals — it only discovers and copies files.
6. **Never touch another module file** unless the user explicitly asked to change more
   than one tab.
7. **Test in isolation**: `useShinyCellModular(data_type = "<type>", enabled_tabs =
   "<tab_id>", ...)`, then `runApp()`.

If the requested feature would require changing the server function's argument list,
prefer adding a new **optional** argument with a default over changing an existing one.

---

## Module File Structure

Every module follows the same layout — **do not deviate from this**:

```r
# Short description of what this tab does
# id     = "<tab_id>"
# title  = "<Human Readable Title>"

############################################### Functions ############################################

# helper / plot functions …
# function names should be unique enough not to collide across modules
# prefer sc_* prefix for shared-looking helpers defined locally

############################################### UI ###################################################

<tab_id>_ui <- function(id, sc1conf, sc1def, ...) {
  ns <- NS(id)
  # shiny UI — use ns() for all input/output IDs
  # sidebarLayout + sidebarPanel / mainPanel is the standard pattern
}

############################################### Server ##############################################

<tab_id>_server <- function(id, sc1conf, sc1def, sc1meta, sc1gene, inpH5, dir_inputs, ...) {
  moduleServer(id, function(input, output, session) {
    # reactive logic, renderPlot, renderDT, downloadHandler …
  })
}

############################################### Registration #########################################

register_tab(
  id          = "<tab_id>",
  title       = "<Human Readable Title>",
  ui          = <tab_id>_ui,
  server      = <tab_id>_server,
  author      = "<your name>",
  description = "…",
  version     = "1.0",
  date        = "Mon YYYY",
  source      = "…",
  contact     = "<your email>"
)
```

### Naming conventions

- Tab id = filename without `.R`, e.g. `bubble_heatmap.R` → `id = "bubble_heatmap"`
- UI function = `<id>_ui`
- Server function = `<id>_server`
- Multi-dataset variants live alongside the single-dataset version: `bubble_heatmap_multi.R`

---

## Server Function Signature

Standard arguments passed to every module server:

| Argument | Type | Contents |
|---|---|---|
| `sc1conf` | `data.table` | Column metadata: `ID`, `UI`, `fCL`, `grp`, `fInt`, `fShow` |
| `sc1def` | named list | Default selections (dimred X/Y, colour, etc.) |
| `sc1meta` | `data.table` | Per-cell metadata; columns match `sc1conf$ID` |
| `sc1gene` | named int vector | Gene name → row index in the H5 file |
| `inpH5` | character | Path to `sc1counts.h5` |
| `dir_inputs` | character | Directory containing all ShinyCell output files |

Additional ATAC-specific args passed when the tab is in `RNA_ATAC`:

| Argument | Contents |
|---|---|
| `sc1conf_atac` | ATAC assay config data.table |
| `sc1meta_atac` | ATAC per-cell metadata |
| `sc1gene_atac` | ATAC peak → row index |
| `inpH5_atac` | Path to ATAC H5 file |

Multi-dataset (`_multi`) tabs receive `sc1conf`, `sc1meta`, etc. as **named lists** keyed by
dataset index (character "1", "2", …).

---

## Files Written by `prepShinyCellModular()`

| File | Contents |
|---|---|
| `sc1conf.rds` | Column config data.table |
| `sc1def.rds` | Default UI selections |
| `sc1meta.fst` or `.rds` | Cell metadata |
| `sc1gene.rds` | Gene → H5 row index |
| `sc1counts.h5` | Sparse count matrix (`grp/data`) |
| `markergenes_lists.parquet` | Marker genes (if `do_markers = TRUE`) |
| `ATAC/sc1conf_atac.rds` etc. | ATAC equivalents (if ATAC assay present) |
| `ATAC/sc1fragmentpaths.rds` | Fragment file paths + cell lists |
| `*_umap3d` reduction | Stored in Seurat, coords saved to metadata (if `do_umap3d = TRUE`) |

---

## Key Patterns in Existing Code

### Reading H5 counts
```r
h5file <- H5File$new(inpH5, mode = "r")
on.exit(try(h5file$close_all(), silent = TRUE), add = TRUE)
h5data <- h5file[["grp"]][["data"]]
vals   <- h5data$read(args = list(sc1gene[gene_name], quote(expr=)))
```

### Subsetting cells
```r
# filter ggData to selected sub-group values only
if (length(inpsub2) != 0 && length(inpsub2) != nlevels(ggData$sub)) {
  ggData <- ggData[sub %in% inpsub2]
}
```

### Resolving a metadata column
```r
# from sc1conf UI label → actual column name in sc1meta
col_id <- sc1conf[UI == inp_label]$ID
```

### sctheme / sList / cList
These globals are defined in the generated `app.R` and are available to all modules at runtime:
- `sctheme(base_size, Xang, XjusH, XYval)` — standard ggplot theme
- `sList` — named vector of font sizes
- `cList` — named list of colour palettes

---

## Adding a New Tab — Checklist

1. Copy the closest existing module as a starting point.
2. Rename the file to `<new_id>.R` and place it under `inst/modules/<data_type>/`.
3. Rename all functions: `<old_id>_ui` → `<new_id>_ui`, etc.
4. Update the `register_tab()` call at the bottom (id, title, description, date).
5. Do **not** touch `useShinyCellModular.R` — it discovers modules automatically.
6. If the tab needs a new file written by prep, add the writing logic to
   `prepShinyCellModular.R` under a clearly labelled section and a `do_*` flag.
7. Test with `useShinyCellModular(data_type = "<type>", enabled_tabs = "<new_id>", …)`.

---

## Adding a `_multi` Variant

Multi variants are structurally identical to single-dataset tabs except:

- Filename: `<id>_multi.R` → automatically picked up for multi-dataset apps.
- `sc1conf`, `sc1meta`, `sc1gene`, `inpH5` are all **named lists**.
- UI typically adds a dataset selector (`selectInput` bound to `names(sc1conf)`).
- Server loops or reacts over the selected dataset index.

---

## Comment Style

Follow the style in the existing files exactly:

```r
# Short sentence comment — no full stop, lowercase after the dash
############################################### Section Header #######################################
# multi-line explanation
# continues here
```

- Section headers use the long `###…###` banner with the label centred (approximately).
- Inline comments are short, lowercase, no trailing punctuation.
- Do **not** add Roxygen `#'` blocks inside module files (those belong in `R/` package files only).

---

## Messages and Warnings

- Use `message(...)` for informational output (not `cat`, not `print`).
- Use `warning(..., call. = FALSE)` for non-fatal issues.
- Use `stop(..., call. = FALSE)` for fatal errors.
- **Do not change the wording of existing messages/warnings** unless explicitly asked.
- New messages should match the terse, direct tone of existing ones, e.g.:
  `message("Copied ", sum(ok), " module(s) into: ", dst_modules)`

---

## Developer Rules (important)

- **Minimum-diff principle**: add only what is asked; do not refactor, rename, or reformat
  working code.
- **No wording changes**: existing `message()`, `warning()`, `stop()`, UI label strings,
  and comments are frozen unless the user asks to change them.
- **No new dependencies** without asking first — the package dependency list is intentionally
  managed in `DESCRIPTION` and the template `app.R`.
- **Preserve existing function signatures** — adding optional arguments with defaults is fine;
  removing or reordering existing arguments is not.
- **Never touch a module file the user didn't ask about.** Each tab is independent by
  design — that isolation is the whole point of the architecture, and an AI agent should
  respect it as strictly as a human contributor would.
- When in doubt about where something belongs, search the existing files first before
  inventing a new pattern.

---

## Troubleshooting Steps

When something is broken, follow this order:

1. **Read the error message literally** — it usually tells you exactly what is missing or wrong; work backwards from there
2. **Identify what changed** — if it worked before, the bug is almost always in the diff
3. **Trace execution order** — does X exist when Y tries to use it; what runs first in the loop
4. **Trace data flow** — where is a variable set, where is it read, does it still hold the right value at that point
5. **Trace logic inside functions** — what does the function assume about its inputs and the environment
6. **Check assumptions** — directory exists, non-null, correct type, correct path
7. **Smallest reproducible case** — isolate the minimum code that triggers the failure

Do not hypothesise about causes before exhausting steps 1 and 2.
