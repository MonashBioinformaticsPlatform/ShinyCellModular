# ShinyCellModular: Developer Guide

This guide is for anyone extending ShinyCellModular: adding a tab, modifying an
existing one, or adding a new data type. It covers the core functions, where
files live, and how to use an AI coding agent safely to speed up tab
development.

---

## 1. The core functions

ShinyCellModular is built around three functions. Two you call yourself;
one you call *inside* a module file.

### `prepShinyCellModular()`

Takes a Seurat object and writes everything a module might need to disk:
config tables, metadata, the HDF5 count matrix, marker genes, 3D reductions,
and (for multiome/ATAC) fragment paths and peak links.

```r
prepShinyCellModular(
  seurat_rds      = "seurat_object.rds",   # or seurat_obj = <object in memory>
  out_dir         = "testing_data_RNA",
  assays_selected = "RNA",
  do_umap3d       = TRUE,   # writes a 3D reduction, needed by cellinfo3D tabs
  do_markers      = TRUE    # writes markergenes_lists.parquet, needed by marker tabs
)
```

If a tab you're adding or extending needs extra files beyond what prep
already writes, drop those files in `out_dir` where the tab can find them
(alongside `sc1conf.rds`, `sc1counts.h5`, etc.), and add the code to read
them in the module's server function. Further changes to
`prepShinyCellModular()` itself can be considered later if you want that
file-generation step automated rather than dropped in manually.

### `useShinyCellModular()`

Generates the actual Shiny app: writes `app.R` and copies the selected
modules into a `modules/` folder next to it.

```r
useShinyCellModular(
  out_dir           = "testing_data/",
  data_type         = "RNA",
  enabled_tabs      = c("cellinfo_cellinfo", "violin_boxplot", "pseudobulk"),
  overwrite_modules = TRUE,   # replaces modules/ entirely, see warning below
  app_title         = "Testing"
)
```

- `data_type` selects which module folder to pull from. This isn't a fixed
  list: `useShinyCellModular()` checks the folder structure under
  `inst/modules/` recursively, so check what folders currently exist there
  rather than assuming the ones shown below. You can add your own
  `data_type` folder simply by creating it.
- `enabled_tabs` is optional; omit it to include every available tab for
  that data type, or list specific tab ids if you want a mix of tabs from
  across several data types.
- **`overwrite_modules = TRUE` replaces the whole `modules/` folder.** If
  you've hand-edited a copied module in place, that edit is gone. Either
  keep `overwrite_modules = FALSE` once you've started customising a
  deployed app, or make your changes in the package source
  (`inst/modules/`) instead of the generated copy, so regenerating the app
  doesn't lose them.

You never need to touch `useShinyCellModular.R` itself to add a tab; it
discovers modules by listing the `data_type` folder at run time. Dropping a
correctly-structured file into that folder is enough for it to show up.

### `register_tab()`

Called at the bottom of every module file. It's how a tab announces itself
to `useShinyCellModular()`, and there is no central list of tabs to edit.

```r
register_tab(
  id          = "violin_boxplot",
  title       = "Violin / BoxPlot",
  ui          = violin_boxplot_ui,
  server      = violin_boxplot_server,
  author      = "Your Name",
  description = "Violin and boxplots for gene expression or metadata",
  version     = "1.0",
  date        = "Jul 2026",
  source      = "internal",
  contact     = "your.email@monash.edu"
)
```

`id` must match the filename (without `.R`): this is what `enabled_tabs`
refers to and how the module is found on disk.

---

## 2. Where module files go

```
inst/modules/
└── <data_type>/   ← e.g. RNA, or any folder you create yourself
```

`useShinyCellModular()` checks this structure recursively rather than
matching against a fixed set of names, so the folders that exist today may
not match what's shown here. Check `inst/modules/` yourself to see current
data types, and create a new folder there if you're adding one.

One file per tab, named `<tab_id>.R`, placed in the folder matching the
data type it applies to. The filename minus `.R` **is** the tab id: it's
what you pass to `enabled_tabs` and what `register_tab(id = ...)` must
match exactly.

Multi-dataset variants live alongside the single-dataset version in the
same folder: `bubble_heatmap.R` and `bubble_heatmap_multi.R` sit side by
side in `RNA/`.

That's the whole placement rule. Nothing else needs to be registered,
imported, or listed anywhere else in the package.

---

## 3. Every module follows the same shape

```r
# Short description of what this tab does
# id     = "<tab_id>"
# title  = "<Human Readable Title>"

############################### Functions ###############################
# helper / plot functions; prefer an sc_* prefix for shared-looking helpers

############################### UI #########################################
<tab_id>_ui <- function(id, sc1conf, sc1def, ...) {
  ns <- NS(id)
  # sidebarLayout + sidebarPanel / mainPanel is the standard pattern
}

############################### Server ######################################
<tab_id>_server <- function(id, sc1conf, sc1def, sc1meta, sc1gene, inpH5, dir_inputs, ...) {
  moduleServer(id, function(input, output, session) {
    # reactive logic, renderPlot, renderDT, downloadHandler …
  })
}

############################### Registration ################################
register_tab(id = "<tab_id>", title = "<Human Readable Title>", ui = <tab_id>_ui,
             server = <tab_id>_server, author = "...", description = "...",
             version = "1.0", date = "Mon YYYY", source = "...", contact = "...")
```

Standard server arguments (all tabs get these; ATAC tabs additionally get
`sc1conf_atac`, `sc1meta_atac`, `sc1gene_atac`, `inpH5_atac`):

| Argument | Contents |
|---|---|
| `sc1conf` | Column metadata table (`ID`, `UI`, `fCL`, `grp`, `fInt`, `fShow`) |
| `sc1def` | Default UI selections |
| `sc1meta` | Per-cell metadata |
| `sc1gene` | Gene name → row index in the H5 file |
| `inpH5` | Path to `sc1counts.h5` |
| `dir_inputs` | Directory containing all prep output |

`sctheme()`, `sList`, and `cList` are defined in the generated `app.R` and
available to every module at runtime for consistent styling.

---

## 4. Using an AI agent to add features to existing tabs

Because every tab is self-contained, an AI coding agent (Claude Code,
Cursor, GitHub Copilot, ChatGPT, or any other agent capable of reading and
editing files) can safely work on one tab without needing to understand
the rest of the app; that's the whole point of the module architecture.
The attached file, **`shinycellmodular-skill/SKILL.md`**, is written
specifically to brief an AI agent on these conventions before it touches
any code.

**The skill file is conventions only, not the codebase.** It tells an
agent the structure, naming rules, and constraints to follow; it does not
give the agent the actual repository. The agent still needs real access to
the code it's meant to change (and, for the minimum-diff rule to mean
anything, the other files it's meant to leave alone).

### How to use it

The mechanics differ by tool, but the requirement is the same everywhere:
the agent needs both the conventions in `SKILL.md` and real access to the
repository files, not one without the other.

- **Agents that work directly on a local checkout** (Claude Code, Cursor,
  and similar tools): place `SKILL.md` wherever that tool looks for
  project-level instructions (for Claude Code, `.claude/skills/` at the
  repo root, or `~/.claude/skills/` for a personal, cross-project copy),
  and run the agent from inside a checkout of the repo. It then has the
  full codebase on disk, so the skill file only adds the conventions
  layer on top; it can genuinely see, and avoid touching, other modules.
- **Chat-based tools with file upload or project knowledge** (a Claude
  Project, Cowork, a custom GPT, etc.): uploading `SKILL.md` alone is not
  enough. Also attach the module file you want changed (and, ideally, the
  `inst/modules/` folder as a whole), or point the agent at the repo
  directly: `https://github.com/MonashBioinformaticsPlatform/ShinyCellModular`.
  Without that, an instruction like "don't touch any other module" is
  meaningless since the agent never had those files in front of it to
  begin with.
- **Any other agent**: paste `SKILL.md` in as a system/context message, and
  make sure the agent has read or fetch access to the actual files it's
  editing, the same way you would for a human contributor working from
  written conventions alone.

### A good request looks like

> "Using the ShinyCellModular conventions in SKILL.md, add a log-scale
> toggle to the attached `violin_boxplot.R` module."

The skill file tells the agent: the exact structure to preserve, which
server arguments are available, the naming conventions, and, critically,
the **minimum-diff rule**: change only what was asked, don't refactor
working code, don't touch `useShinyCellModular.R`, and don't rename
existing functions or arguments. That rule only holds if the agent
actually has the surrounding codebase available to potentially touch; if
it was only ever given one file, there's nothing to "not touch" and the
constraint is doing no real work.

### What to check before merging an AI-assisted change

1. Did it touch only the one module file (or `prepShinyCellModular.R` if
   the feature needs a new prep flag)?
2. Does `register_tab()` still match the filename and existing `id`?
3. Did it preserve existing message/warning wording, or only add new
   messages in the same terse style?
4. Test it in isolation: `useShinyCellModular(data_type = "<type>",
   enabled_tabs = "<tab_id>", ...)`.

---

## 5. Adding a brand-new tab: checklist

1. Copy the closest existing module as a starting point.
2. Rename the file to `<new_id>.R` and place it under
   `inst/modules/<data_type>/`.
3. Rename all functions: `<old_id>_ui` → `<new_id>_ui`, etc.
4. Update the `register_tab()` call (id, title, description, date).
5. Leave `useShinyCellModular.R` alone; it discovers modules automatically.
6. If the tab needs a new file written by prep, add the writing logic to
   `prepShinyCellModular.R` under a clearly labelled section and a `do_*`
   flag, so users can opt in without slowing down prep for tabs they don't
   use.
7. Test with `useShinyCellModular(data_type = "<type>", enabled_tabs =
   "<new_id>", ...)`.
