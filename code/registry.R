# this registry has all possible tabs in the modules folder
tab_registry <- list(
  umap = list(
    title  = "UMAP",
    ui     = mod_umap_ui,
    server = mod_umap_server
  ),
  viobox = list(
    title  = "Violin/Box",
    ui     = mod_viobox_ui,
    server = mod_viobox_server
  ),
  dotplot = list(
    title  = "Dot plot",
    ui     = mod_dotplot_ui,
    server = mod_dotplot_server
  )
)