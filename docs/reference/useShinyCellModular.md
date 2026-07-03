# Generate a modular ShinyCellModular Shiny app

Takes the output files from
[`prepShinyCellModular`](https://monashbioinformaticsplatform.github.io/ShinyCellModular/reference/prepShinyCellModular.md)
and writes a ready-to-run `app.R` with the selected module tabs into
`out_dir`.

## Usage

``` r
useShinyCellModular(
  out_dir,
  shinycellmodular.dir.src = NULL,
  rsconnect.deploy = FALSE,
  data_type = NULL,
  enabled_tabs = NULL,
  overwrite_modules = FALSE,
  disable_ui_server = TRUE,
  app_title = NULL,
  navbar_color = "#007BA7"
)
```

## Arguments

- out_dir:

  Directory containing prepared prepShinyCellModular output files
  (`sc1conf.rds`, `sc1meta.rds`, etc.).

- shinycellmodular.dir.src:

  Path to the ShinyCellModular source directory containing `modules/`.
  Default: `system.file('', package = 'ShinyCellModular')`.

- rsconnect.deploy:

  Write an rsconnect manifest for deployment. Default: `FALSE`.

- data_type:

  Preset tab selection: `'RNA'`, `'RNA_ATAC'`, or `'SPATIAL'`. Default:
  `'RNA'`.

- enabled_tabs:

  Character vector of tab IDs to include. Overrides `data_type` presets.
  Default: `NULL` (uses all tabs for `data_type`).

- overwrite_modules:

  Remove and replace existing `modules/` folder. Default: `FALSE`.

- disable_ui_server:

  Rename legacy `ui.R` and `server.R` to `.bak`. Default: `TRUE`.

- app_title:

  Title shown in the app navbar. Required.

## Value

Invisibly returns `NULL`. Writes `app.R` and `modules/` to `out_dir`.

## Examples

``` r
if (FALSE) { # \dontrun{
useShinyCellModular(
  out_dir  = "my_app_files/",
  data_type  = "RNA",
  app_title  = "My scRNA-seq App"
)
} # }
```
