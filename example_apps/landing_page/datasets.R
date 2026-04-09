# ==============================================================================
# datasets.R
# Configuration for the ShinyCellModular landing page.
# Add / edit entries here to update the portal — do not touch app.R.
#
# Each entry accepts:
#   id          - unique identifier (no spaces)
#   title       - display name of the dataset
#   description - one or two sentence description shown on the card
#   data_type   - one of: "RNA", "ATAC", "Multiome", "Spatial"
#   organism    - e.g. "Human", "Mouse"
#   platform    - e.g. "10x Genomics", "Parse Biosciences"
#   n_cells     - character string, e.g. "14,200"
#   n_clusters  - integer
#   extra_meta  - named character vector of any extra key/value pairs shown on
#                 the card (optional, set to NULL to skip)
#   app_url     - URL to the deployed ShinyCellModular Shiny app
#   rmd_url     - URL to the RMarkdown tutorial for recreating this app
#   updated     - character string shown as last-updated label
# ==============================================================================

datasets <- list(

  # --- scRNA-seq ---------------------------------------------------------------
  list(
    id          = "scrna_pbmc",
    title       = "PBMC scRNA-seq",
    description = "Peripheral blood mononuclear cells profiled with 10x Genomics. Reference dataset for immune cell type annotation and differential abundance.",
    data_type   = "RNA",
    organism    = "Human",
    platform    = "10x Genomics",
    n_cells     = "2,700",
    n_clusters  = 14L,
    extra_meta  = c("genes" = "33k"),
    app_url     = "https://your-server/content/scrna-pbmc",       # TODO: replace URL
    rmd_url     = "https://your-docs/scrna-pbmc-tutorial",        # TODO: replace URL
    updated     = "2 weeks ago"
  ),

  # --- ATAC --------------------------------------------------------------------
  list(
    id          = "atac_breast",
    title       = "Breast cancer ATAC",
    description = "Chromatin accessibility landscape of primary breast tumours across subtypes. Includes peak-level annotation and motif enrichment.",
    data_type   = "ATAC",
    organism    = "Human",
    platform    = "10x ATAC",
    n_cells     = "6,800",
    n_clusters  = 11L,
    extra_meta  = NULL,
    app_url     = "https://your-server/content/atac-breast",      # TODO: replace URL
    rmd_url     = "https://your-docs/atac-breast-tutorial",       # TODO: replace URL
    updated     = "1 month ago"
  ),

  # --- Multiome ----------------------------------------------------------------
  list(
    id          = "multiome_cortex",
    title       = "Mouse cortex multiome",
    description = "Joint ATAC + RNA profiling of mouse primary motor cortex. Chromatin accessibility and gene expression explored side by side.",
    data_type   = "Multiome",
    organism    = "Mouse",
    platform    = "10x Multiome",
    n_cells     = "9,400",
    n_clusters  = 21L,
    extra_meta  = c("modalities" = "RNA + ATAC"),
    app_url     = "https://your-server/content/multiome-cortex",  # TODO: replace URL
    rmd_url     = "https://your-docs/multiome-cortex-tutorial",   # TODO: replace URL
    updated     = "3 weeks ago"
  ),

  # --- Spatial -----------------------------------------------------------------
  list(
    id          = "spatial_liver",
    title       = "Human liver spatial",
    description = "Visium spatial transcriptomics of healthy adult human liver. Spot-level gene expression overlaid on tissue histology.",
    data_type   = "Spatial",
    organism    = "Human",
    platform    = "10x Visium",
    n_cells     = "3,100",
    n_clusters  = 8L,
    extra_meta  = c("sections" = "4 donors"),
    app_url     = "https://your-server/content/spatial-liver",    # TODO: replace URL
    rmd_url     = "https://your-docs/spatial-liver-tutorial",     # TODO: replace URL
    updated     = "5 days ago"
  )

)
