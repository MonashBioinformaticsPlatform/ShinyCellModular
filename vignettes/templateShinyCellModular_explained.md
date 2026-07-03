# templateShinyCellModular()

`templateShinyCellModular()` scaffolds a new module/tab file. It doesn't do any analysis — it writes a skeleton `.R` file following the standard module structure (Functions / UI / Server / Registration) into the right `modules/<data_type>/` subfolder, so a new tab can be built by filling in placeholders instead of copy-pasting an existing module by hand.

---

## Function signature

```r
templateShinyCellModular(
    id,
    title,
    data_type                = "RNA",
    multi                     = FALSE,
    description               = "TODO: describe this tab",
    author                    = NULL,
    contact                   = NULL,
    shinycellmodular.dir.src  = NULL,
    overwrite                 = FALSE
)
```

---

## Arguments

| Argument | Default | Description |
|---|---|---|
| `id` | — (required) | Tab id. Becomes the filename (sans extension) and the prefix for the generated function names (`<id>_ui`, `<id>_server`). Must be unique across `modules/`. |
| `title` | — (required) | Human readable title shown on the tab. |
| `data_type` | `"RNA"` | Which `modules/` subfolder to write into: `"RNA"`, `"RNA_ATAC"`, or `"SPATIAL"`. Validated with `match.arg()`. |
| `multi` | `FALSE` | Write a `_multi` variant. `_multi` is appended to `id` for the filename, function names, and registered tab id, matching the existing multi-dataset naming convention. |
| `description` | `"TODO: describe this tab"` | Short description used in `register_tab()`. |
| `author` | `NULL` | Author name used in `register_tab()`. Required. |
| `contact` | `NULL` | Contact email used in `register_tab()`. Required. |
| `shinycellmodular.dir.src` | `NULL` | Path to the ShinyCellModular source directory containing `modules/`. Auto-detected via `find.package("ShinyCellModular")` if not supplied. |
| `overwrite` | `FALSE` | Overwrite an existing file with the same name. |

**Return value:** invisibly returns the path to the new file. The real output is the written `.R` file under `modules/<data_type>/`.

---

## Step by step

### 1. Required-argument checks
`id`, `title`, `author`, and `contact` are each checked with `missing()` / `is.null()` / `!nzchar()` and stop with a usage example if absent — no defaults are silently substituted for these.

### 2. Validate `data_type`
`match.arg()` restricts `data_type` to `"RNA"`, `"RNA_ATAC"`, or `"SPATIAL"`.

### 3. Resolve `shinycellmodular.dir.src`
Same auto-detect pattern as `useShinyCellModular()`: if not supplied, tries `find.package("ShinyCellModular")` and checks it has a `modules/` folder, else stops with instructions to pass the path explicitly. Prints a message when auto-detected, then normalises the path with `normalizePath(mustWork = TRUE)`.

### 4. Build the target id, folder, and file path
`module_id` is `id` with `_multi` appended when `multi = TRUE` (e.g. `bubble_heatmap` → `bubble_heatmap_multi`), matching the existing naming convention. The destination is `shinycellmodular.dir.src/modules/<data_type>/<module_id>.R`.

- Stops if the target `modules/<data_type>/` folder doesn't exist.
- Stops if a file already exists at that path and `overwrite = FALSE`, telling you to set `overwrite = TRUE` to replace it.

### 5. Fill in the template and write the file
A fixed template string containing the four standard sections (Functions / UI / Server / Registration) has its placeholders substituted with `gsub(..., fixed = TRUE)`:

| Placeholder | Replaced with |
|---|---|
| `__ID__` | `module_id` |
| `__TITLE__` | `title` |
| `__AUTHOR__` | `author` |
| `__CONTACT__` | `contact` |
| `__DESCRIPTION__` | `description` |
| `__DATE__` | `format(Sys.Date(), "%b %Y")` |

The generated file already has correctly namespaced `<id>_ui` / `<id>_server` functions and a `register_tab()` call at the bottom with `version = "1.0"` and `source = "MGBP custom"` pre-filled, so the only work left for the developer is filling in the actual UI/server logic between the section banners.

Writes the file with `writeLines()` and prints a confirmation message with the destination path.

---

## Notes

- Because the template is a single fixed string, any future changes to the standard module layout (e.g. a new required section) need to be made in this template too, or new modules will drift from the convention documented in the module-authoring guide.
- `overwrite = TRUE` replaces the file outright — there's no `.bak` step here (unlike `useShinyCellModular()`'s handling of `ui.R`/`server.R`), so uncommitted work in an existing module file will be lost.
