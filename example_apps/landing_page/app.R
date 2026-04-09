# ==============================================================================
# app.R
# ShinyCellModular landing page — showcases available example apps.
#
# To add or edit datasets edit datasets.R only.
# To change the appearance edit the CSS block inside landing_page_ui() below.
# ==============================================================================

library(shiny)
library(bslib)

source("datasets.R")

# ------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------

# badge colour per data type — Monash secondary palette tints
badge_colours <- list(
  RNA      = list(bg = "#daeef8", text = "#027EB6"),   # Monash light blue
  ATAC     = list(bg = "#ede9f5", text = "#4a4680"),   # Monash purple tint
  Multiome = list(bg = "#d6eedf", text = "#005c19"),   # Monash green tint
  Spatial  = list(bg = "#fae3d9", text = "#a02900")    # Monash orange tint
)

make_badge <- function(data_type) {
  col <- badge_colours[[data_type]]
  if (is.null(col)) col <- list(bg = "#F1EFE8", text = "#444441")
  tags$span(
    data_type,
    style = paste0(
      "font-size:11px; padding:3px 9px; border-radius:6px; white-space:nowrap;",
      "background:", col$bg, "; color:", col$text, ";"
    )
  )
}

make_meta_pill <- function(label, value) {
  tags$span(
    paste0(label, ": ", value),
    style = "font-size:12px; color:#505050; background:#F6F6F6;
             padding:2px 8px; border-radius:4px;"
  )
}

# build one dataset card
make_card <- function(ds) {
  # extra meta pills (cells + clusters always shown, plus any extras)
  meta_pills <- tagList(
    make_meta_pill("cells",    ds$n_cells),
    make_meta_pill("clusters", ds$n_clusters)
  )
  if (!is.null(ds$extra_meta)) {
    for (nm in names(ds$extra_meta)) {
      meta_pills <- tagAppendChild(meta_pills, make_meta_pill(nm, ds$extra_meta[[nm]]))
    }
  }
  
  div(
    class = "scm-card",
    # --- header ---------------------------------------------------------------
    div(
      style = "display:flex; justify-content:space-between; align-items:flex-start;
               margin-bottom:10px;",
      div(
        tags$p(ds$title,    class = "card-title"),
        tags$p(paste0(ds$platform, " \u00b7 ", ds$organism), class = "card-subtitle")
      ),
      make_badge(ds$data_type)
    ),
    # --- description ----------------------------------------------------------
    tags$p(ds$description, class = "card-desc"),
    # --- meta pills -----------------------------------------------------------
    div(style = "display:flex; flex-wrap:wrap; gap:6px; margin:8px 0;", meta_pills),
    # --- footer ---------------------------------------------------------------
    div(
      class = "card-footer",
      tags$span(paste0("updated ", ds$updated), class = "card-updated"),
      div(
        style = "display:flex; gap:8px;",
        tags$a(
          href   = ds$rmd_url,
          target = "_blank",
          class  = "scm-btn scm-btn-secondary",
          "How-to guide"
        ),
        tags$a(
          href   = ds$app_url,
          target = "_blank",
          class  = "scm-btn scm-btn-primary",
          "Open app \u2197"
        )
      )
    )
  )
}

# ------------------------------------------------------------------------------
# ui
# ------------------------------------------------------------------------------

ui <- fluidPage(
  
  # page title for the browser tab
  title = "ShinyCellModular Portal",
  
  tags$head(tags$style(HTML("

    /* ---- reset & base ---------------------------------------------------- */
    /* Monash primary palette: blue #006DAE, black #000000, dark grey #3c3c3c,
       grey #505050, light grey #F6F6F6                                        */
    * { box-sizing: border-box; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
      background: #F6F6F6;
      color: #000000;
      margin: 0;
    }

    /* ---- top nav bar — Monash blue ---------------------------------------- */
    .monash-nav {
      background: #006DAE;
      padding: 0 1.5rem;
      display: flex; align-items: center; gap: 12px;
      height: 52px;
    }
    .monash-nav-wordmark {
      font-size: 15px; font-weight: 700; color: #fff;
      letter-spacing: 0.08em; text-transform: uppercase;
    }
    .monash-nav-divider {
      width: 1px; height: 22px; background: rgba(255,255,255,0.35);
    }
    .monash-nav-sub {
      font-size: 12px; color: rgba(255,255,255,0.85); letter-spacing: 0.03em;
    }

    /* ---- bootstrap reset — prevent row/col interfering with our grid ------ */
    .portal-wrap .row { margin-left: 0; margin-right: 0; }
    .portal-wrap .col, .portal-wrap [class*='col-'] {
      padding-left: 0; padding-right: 0; float: none;
    }

    /* ---- page layout ----------------------------------------------------- */
    .portal-wrap   { max-width: 960px; margin: 0 auto; padding: 2.5rem 1.5rem; }

    /* ---- header ---------------------------------------------------------- */
    .portal-header { margin-bottom: 2rem; }
    .portal-title  {
      font-size: 24px; font-weight: 700; color: #000000; margin: 0 0 4px;
    }
    /* blue underline accent on the title */
    .portal-title span {
      border-bottom: 3px solid #006DAE; padding-bottom: 2px;
    }
    .portal-sub    { font-size: 14px; color: #505050; margin: 0; }

    /* ---- controls row ----------------------------------------------------- */
    .controls-row  {
      display: flex; flex-wrap: wrap; gap: 10px;
      align-items: center; margin-bottom: 1.5rem;
    }
    .controls-row input[type=text] {
      padding: 7px 12px; font-size: 13px;
      border: 1px solid #c8c8c8; border-radius: 4px;
      background: #fff; color: #000000; width: 220px;
      outline: none;
    }
    .controls-row input[type=text]:focus { border-color: #006DAE; }

    /* ---- filter tags ------------------------------------------------------ */
    .tag-group    { display: flex; gap: 6px; flex-wrap: wrap; }
    .filter-tag   {
      font-size: 12px; padding: 5px 12px;
      border: 1px solid #c8c8c8; border-radius: 4px;
      background: #fff; color: #505050;
      cursor: pointer; user-select: none; transition: all 0.1s;
    }
    .filter-tag:hover   { border-color: #006DAE; color: #006DAE; }
    .filter-tag.active  {
      background: #006DAE; border-color: #006DAE; color: #fff;
    }

    /* ---- stats strip ------------------------------------------------------ */
    .stats-strip  {
      display: grid; grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 12px; margin-bottom: 2rem;
    }
    .stat-box     {
      background: #fff; border-top: 3px solid #006DAE;
      border-left: 1px solid #e0e0e0;
      border-right: 1px solid #e0e0e0;
      border-bottom: 1px solid #e0e0e0;
      border-radius: 0 0 4px 4px; padding: 14px 16px;
    }
    .stat-label   { font-size: 12px; color: #505050; text-transform: uppercase;
                    letter-spacing: 0.04em; }
    .stat-value   { font-size: 24px; font-weight: 700; color: #006DAE; margin-top: 2px; }

    /* ---- section label ---------------------------------------------------- */
    .section-label {
      font-size: 11px; font-weight: 700; color: #505050;
      letter-spacing: 0.08em; text-transform: uppercase;
      margin-bottom: 14px;
    }

    /* ---- card grid -------------------------------------------------------- */
    /* layout is handled server-side with fluidRow + column(6)               */
    .card-grid    { width: 100%; }
    .card-grid .row { margin-left: -8px; margin-right: -8px; }
    .card-grid .row > div { padding-left: 8px; padding-right: 8px;
                            margin-bottom: 16px; }
    .scm-card     {
      background: #fff;
      border: 1px solid #e0e0e0;
      border-top: 3px solid #3c3c3c;
      border-radius: 0 0 4px 4px; padding: 1.25rem;
      display: flex; flex-direction: column;
      transition: border-top-color 0.15s;
    }
    .scm-card:hover { border-top-color: #006DAE; }

    /* card typography */
    .card-title    { font-size: 15px; font-weight: 700; color: #000000; margin: 0 0 2px; }
    .card-subtitle { font-size: 12px; color: #505050; margin: 0; }
    .card-desc     { font-size: 13px; color: #3c3c3c; line-height: 1.55; margin: 0; }

    /* card footer */
    .card-footer   {
      margin-top: auto; padding-top: 12px;
      border-top: 1px solid #e0e0e0;
      display: flex; justify-content: space-between; align-items: center;
    }
    .card-updated  { font-size: 11px; color: #969696; }

    /* buttons */
    .scm-btn       {
      font-size: 12px; padding: 6px 13px;
      border-radius: 4px; text-decoration: none;
      cursor: pointer; white-space: nowrap;
      transition: background 0.1s, color 0.1s;
    }
    /* primary — Monash blue fill */
    .scm-btn-primary {
      background: #006DAE; color: #fff; border: 1px solid #006DAE;
    }
    .scm-btn-primary:hover  { background: #005a8e; color: #fff; border-color: #005a8e; }
    /* secondary — outlined */
    .scm-btn-secondary {
      background: transparent; color: #006DAE;
      border: 1px solid #006DAE;
    }
    .scm-btn-secondary:hover { background: #daeef8; color: #006DAE; }

    /* hidden cards when filtered out */
    .scm-card.hidden { display: none; }

    /* no-results message */
    #no-results {
      display: none; grid-column: 1 / -1;
      font-size: 14px; color: #505050; padding: 2rem 0;
    }

  "))),
  
  # --- Monash nav bar -----------------------------------------------------------
  div(
    class = "monash-nav",
    tags$span("Monash", class = "monash-nav-wordmark"),
    div(class = "monash-nav-divider"),
    tags$span("Genomics and Bioinformatics Platform", class = "monash-nav-sub")
  ),
  
  div(
    class = "portal-wrap",
    
    # --- header ---------------------------------------------------------------
    div(
      class = "portal-header",
      tags$h1(tags$span("ShinyCellModular"), " Portal", class = "portal-title"),
      tags$p(
        "Example apps by data type — select a dataset to explore or follow the how-to guide to build your own",
        class = "portal-sub"
      )
    ),
    
    # --- controls row ---------------------------------------------------------
    div(
      class = "controls-row",
      tags$input(
        type        = "text",
        id          = "search-box",
        placeholder = "Search datasets\u2026",
        oninput     = "filterCards()"
      ),
      div(
        class = "tag-group",
        id    = "tag-group",
        # "All" tag starts active; others generated from unique data_types
        tags$div("All",      class = "filter-tag active", `data-type` = "All",
                 onclick = "setTag(this)"),
        tags$div("RNA",      class = "filter-tag", `data-type` = "RNA",
                 onclick = "setTag(this)"),
        tags$div("ATAC",     class = "filter-tag", `data-type` = "ATAC",
                 onclick = "setTag(this)"),
        tags$div("Multiome", class = "filter-tag", `data-type` = "Multiome",
                 onclick = "setTag(this)"),
        tags$div("Spatial",  class = "filter-tag", `data-type` = "Spatial",
                 onclick = "setTag(this)")
      )
    ),
    
    # --- stats strip ----------------------------------------------------------
    # built statically — values derived from datasets at startup
    local({
      n_datasets  <- length(datasets)
      total_cells <- sum(sapply(datasets, function(d) as.integer(gsub(",", "", d$n_cells))))
      data_types  <- length(unique(sapply(datasets, `[[`, "data_type")))
      div(
        style = "display:flex; gap:12px; margin-bottom:2rem;",
        div(class = "stat-box", style = "flex:1;",
            tags$p("datasets",   class = "stat-label"),
            tags$p(n_datasets,   class = "stat-value")),
        div(class = "stat-box", style = "flex:1;",
            tags$p("total cells", class = "stat-label"),
            tags$p(formatC(total_cells, format = "d", big.mark = ","), class = "stat-value")),
        div(class = "stat-box", style = "flex:1;",
            tags$p("data types", class = "stat-label"),
            tags$p(data_types,   class = "stat-value"))
      )
    }),
    
    # --- card grid ------------------------------------------------------------
    # cards built statically at startup — no uiOutput wrapper to break the grid
    tags$p("available datasets", class = "section-label"),
    div(
      id    = "card-grid",
      style = "display:flex; flex-wrap:wrap; gap:16px;",
      tagList(lapply(datasets, function(ds) {
        card <- make_card(ds)
        card$attribs[["data-type"]]  <- ds$data_type
        card$attribs[["style"]]      <- "width:calc(50% - 8px); min-width:220px;"
        card
      })),
      tags$div(id = "no-results", style = "display:none; width:100%;",
               "No datasets match your search.")
    )
  ),
  
  # --- client-side filter logic (no server round-trip needed) -----------------
  tags$script(HTML("

    var activeTag = 'All';

    function setTag(el) {
      document.querySelectorAll('.filter-tag').forEach(function(t) {
        t.classList.remove('active');
      });
      el.classList.add('active');
      activeTag = el.getAttribute('data-type');
      filterCards();
    }

    function filterCards() {
      var query = document.getElementById('search-box').value.toLowerCase();
      var cards = document.querySelectorAll('.scm-card');
      var visible = 0;

      cards.forEach(function(card) {
        var typeMatch = (activeTag === 'All') ||
                        (card.getAttribute('data-type') === activeTag);
        var textMatch = card.textContent.toLowerCase().indexOf(query) !== -1;

        if (typeMatch && textMatch) {
          card.classList.remove('hidden');
          visible++;
        } else {
          card.classList.add('hidden');
        }
      });

      var noRes = document.getElementById('no-results');
      noRes.style.display = (visible === 0) ? 'block' : 'none';
    }

  "))
)

# ------------------------------------------------------------------------------
# server
# ------------------------------------------------------------------------------

server <- function(input, output, session) {
  # stats strip and cards are rendered statically in ui — nothing to do here
}


shinyApp(ui, server)