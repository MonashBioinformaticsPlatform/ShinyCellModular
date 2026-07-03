# Create a template file for a new ShinyCellModular module/tab

Writes a skeleton .R file following the standard module structure
(Functions / UI / Server / Registration) used across ShinyCellModular
modules, so a new tab can be built by filling in the placeholders
instead of copy-pasting an existing module by hand.

## Usage

``` r
templateShinyCellModular(
  id,
  title,
  data_type = "RNA",
  multi = FALSE,
  description = "TODO: describe this tab",
  author = NULL,
  contact = NULL,
  shinycellmodular.dir.src = NULL,
  overwrite = FALSE
)
```

## Arguments

- id:

  Tab id — used as the filename (sans extension) and as the prefix for
  the generated function names (`<id>_ui`, `<id>_server`). Must be
  unique across `modules/`. Required.

- title:

  Human readable title shown on the tab. Required.

- data_type:

  Which `modules/` subfolder to write into: `'RNA'`, `'RNA_ATAC'` or
  `'SPATIAL'`. Default: `'RNA'`.

- multi:

  Write a `_multi` variant. `_multi` is appended to `id` for the
  filename, function names and registered tab id, per the existing
  multi-dataset naming convention. Default: `FALSE`.

- description:

  Short description used in `register_tab()`. Default:
  `"TODO: describe this tab"`.

- author:

  Author name used in `register_tab()`. Required.

- contact:

  Contact email used in `register_tab()`. Required.

- shinycellmodular.dir.src:

  Path to the ShinyCellModular source directory containing `modules/`.
  Default: `system.file('', package = 'ShinyCellModular')`.

- overwrite:

  Overwrite an existing file with the same name. Default: `FALSE`.

## Value

Invisibly returns the path to the new file.

## Examples

``` r
if (FALSE) { # \dontrun{
templateShinyCellModular(
  id      = "my_new_tab",
  title   = "My New Tab",
  data_type = "RNA",
  author  = "Your Name",
  contact = "your.email@monash.edu"
)
} # }
```
