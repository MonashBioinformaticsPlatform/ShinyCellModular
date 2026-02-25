# id     = "pseudobulk"
# title  = "Pseudobulk DE (voom limma)"

############################################### Functions ###########################################

sc_has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

sc_resolve_meta_col <- function(sc1conf, sc1meta, x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(NULL)
  
  if (!is.null(sc1conf) && "UI" %in% names(sc1conf) && x %in% sc1conf$UI) {
    id <- sc1conf[UI == x]$ID
    if (length(id) == 1 && id %in% names(sc1meta)) return(id)
  }
  
  if (x %in% names(sc1meta)) return(x)
  
  NULL
}

sc_h5_find_dataset <- function(h5, candidates) {
  for (p in candidates) {
    ok <- FALSE
    obj <- NULL
    ok <- tryCatch({
      obj <- h5[[p]]
      !is.null(obj)
    }, error = function(e) FALSE)
    if (isTRUE(ok)) return(list(path = p, obj = obj))
  }
  NULL
}

sc_h5_read_counts_block <- function(h5_path, gene_idx, cell_idx,
                                    dataset_candidates = c("counts/counts", "counts/data", "grp/data")) {
  shiny::validate(
    shiny::need(sc_has_pkg("hdf5r"), "Package hdf5r is required to read H5 counts.")
  )
  
  h5 <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  
  found <- sc_h5_find_dataset(h5, dataset_candidates)
  shiny::validate(
    shiny::need(!is.null(found), paste0(
      "Could not find a counts dataset in ", basename(h5_path), ". ",
      "Tried: ", paste(dataset_candidates, collapse = ", ")
    ))
  )
  
  dset <- found$obj
  dset[gene_idx, cell_idx, drop = FALSE]
}

sc_make_pseudobulk <- function(sc1conf,
                               sc1meta,
                               sc1gene,
                               h5_counts_path,
                               keep_cells_idx,
                               unit_col_ui_or_meta,
                               condition_col_ui_or_meta,
                               covar_cols_ui_or_meta = character(0),
                               rep_within_cond = FALSE,
                               chunk_genes = 1000,
                               dataset_candidates = c("counts/counts", "counts/data", "grp/data")) {
  
  unit_col <- sc_resolve_meta_col(sc1conf, sc1meta, unit_col_ui_or_meta)
  cond_col <- sc_resolve_meta_col(sc1conf, sc1meta, condition_col_ui_or_meta)
  
  shiny::validate(shiny::need(!is.null(unit_col), paste0("Cannot resolve replicate column: ", unit_col_ui_or_meta)))
  shiny::validate(shiny::need(!is.null(cond_col), paste0("Cannot resolve condition column: ", condition_col_ui_or_meta)))
  shiny::validate(shiny::need(length(keep_cells_idx) > 0, "No cells left after filtering."))
  
  unit_vec <- as.character(sc1meta[[unit_col]][keep_cells_idx])
  cond_vec <- as.character(sc1meta[[cond_col]][keep_cells_idx])
  
  ok <- !is.na(unit_vec) & nzchar(unit_vec) & !is.na(cond_vec) & nzchar(cond_vec)
  keep_cells_idx <- keep_cells_idx[ok]
  unit_vec <- unit_vec[ok]
  cond_vec <- cond_vec[ok]
  
  shiny::validate(shiny::need(length(keep_cells_idx) > 0, "No cells left after removing NA replicate or condition."))
  
  unit_key <- unit_vec
  if (isTRUE(rep_within_cond)) {
    unit_key <- paste0(unit_vec, "___", cond_vec)
  }
  
  uniq_units <- unique(unit_key)
  shiny::validate(shiny::need(length(uniq_units) > 1, "Only 1 pseudobulk unit present after filtering."))
  
  unit_index <- match(unit_key, uniq_units)
  
  unit_condition <- tapply(cond_vec, unit_key, function(x) unique(x[!is.na(x) & nzchar(x)]))
  ok_one <- vapply(unit_condition, function(x) length(x) == 1, logical(1))
  shiny::validate(
    shiny::need(all(ok_one),
                "Some replicates map to multiple conditions. Check replicate and condition columns, or enable replicate within condition.")
  )
  
  unit_condition <- vapply(unit_condition, function(x) x[[1]], character(1))
  unit_condition <- unit_condition[uniq_units]
  
  covar_cols <- vapply(covar_cols_ui_or_meta, function(z) sc_resolve_meta_col(sc1conf, sc1meta, z), character(1))
  covar_cols <- covar_cols[!is.na(covar_cols) & nzchar(covar_cols)]
  covar_cols <- unique(covar_cols)
  
  covar_df <- NULL
  if (length(covar_cols) > 0) {
    covar_df <- lapply(covar_cols, function(col) {
      v <- sc1meta[[col]][keep_cells_idx]
      tapply(v, unit_key, function(x) {
        ux <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
        if (length(ux) == 0) NA_character_
        else if (length(ux) == 1) ux[[1]]
        else NA_character_
      })[uniq_units]
    })
    covar_df <- as.data.frame(covar_df, stringsAsFactors = FALSE)
    names(covar_df) <- covar_cols
    
    bad <- vapply(covar_df, function(x) any(is.na(x)), logical(1))
    shiny::validate(
      shiny::need(!any(bad), paste0("Some covariates are not constant within unit: ", paste(names(covar_df)[bad], collapse = ", ")))
    )
  }
  
  gene_idx <- as.integer(unname(sc1gene))
  gene_names <- names(sc1gene)
  shiny::validate(shiny::need(length(gene_idx) > 1, "sc1gene is empty or invalid."))
  
  nG <- length(gene_idx)
  nU <- length(uniq_units)
  
  pb <- matrix(0, nrow = nG, ncol = nU)
  rownames(pb) <- gene_names
  colnames(pb) <- uniq_units
  
  gene_chunks <- split(seq_len(nG), ceiling(seq_len(nG) / chunk_genes))
  
  do_block <- function(ii) {
    gi <- gene_idx[ii]
    m <- sc_h5_read_counts_block(
      h5_counts_path,
      gene_idx = gi,
      cell_idx = keep_cells_idx,
      dataset_candidates = dataset_candidates
    )
    
    if (!is.matrix(m)) m <- as.matrix(m)
    if (nrow(m) != length(ii)) {
      m <- matrix(m, nrow = length(ii))
    }
    
    rs <- rowsum(t(m), group = unit_index, reorder = FALSE)
    pb[ii, ] <<- t(rs)
    NULL
  }
  
  shiny::withProgress(message = "Building pseudobulk counts", value = 0, {
    nC <- length(gene_chunks)
    for (k in seq_len(nC)) {
      shiny::incProgress(1 / nC, detail = paste0("Chunk ", k, " of ", nC))
      do_block(gene_chunks[[k]])
    }
  })
  
  design_df <- data.frame(
    unit = uniq_units,
    condition = unit_condition,
    stringsAsFactors = FALSE
  )
  if (!is.null(covar_df) && ncol(covar_df) > 0) {
    design_df <- cbind(design_df, covar_df)
  }
  
  list(
    pseudobulk = pb,
    design_df = design_df,
    keep_cells_idx = keep_cells_idx,
    unit_ids = uniq_units
  )
}

sc_run_voom_limma <- function(pseudobulk,
                              design_df,
                              model_str = "~ condition",
                              contrast_str = "",
                              robust_ebayes = TRUE,
                              reference_level = NULL) {
  
  shiny::validate(shiny::need(sc_has_pkg("edgeR"), "edgeR is required."))
  shiny::validate(shiny::need(sc_has_pkg("limma"), "limma is required."))
  
  design_df <- as.data.frame(design_df, stringsAsFactors = FALSE)
  shiny::validate(shiny::need("condition" %in% names(design_df), "design_df must contain condition."))
  
  design_df$condition <- factor(design_df$condition)
  shiny::validate(shiny::need(nlevels(design_df$condition) > 1, "Condition has only 1 level after filtering."))
  
  if (!is.null(reference_level) && nzchar(reference_level)) {
    shiny::validate(shiny::need(reference_level %in% levels(design_df$condition), "Selected reference is not in condition levels."))
    design_df$condition <- stats::relevel(design_df$condition, ref = reference_level)
  }
  
  fml <- stats::as.formula(model_str)
  
  dge <- edgeR::DGEList(pseudobulk)
  dge <- edgeR::calcNormFactors(dge)
  
  design <- stats::model.matrix(fml, data = design_df)
  
  v <- limma::voom(dge, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design = design)
  
  if (!is.null(contrast_str) && nzchar(contrast_str)) {
    contr <- limma::makeContrasts(contrasts = contrast_str, levels = design)
    fit <- limma::contrasts.fit(fit, contr)
  }
  
  fit <- limma::eBayes(fit, robust = isTRUE(robust_ebayes))
  
  tt <- limma::topTable(fit, n = Inf, adjust.method = "BH", sort.by = "P")
  tt$gene <- rownames(tt)
  
  list(
    dge = dge,
    voom = v,
    design = design,
    fit = fit,
    table = tt,
    design_df = design_df
  )
}


sc_pseudobulk_plots <- function(tt, p_cut = 0.05, fc_up = 1, fc_down = -1) {

  shiny::validate(shiny::need(sc_has_pkg("ggplot2"), "ggplot2 is required for plotting."))

  tt <- as.data.frame(tt, stringsAsFactors = FALSE)

  if (!("gene" %in% names(tt))) tt$gene <- rownames(tt)
  if (!("AveExpr" %in% names(tt))) tt$AveExpr <- NA_real_
  if (!("logFC" %in% names(tt))) tt$logFC <- NA_real_
  if (!("adj.P.Val" %in% names(tt))) tt$adj.P.Val <- NA_real_

  ggData <- data.table::as.data.table(tt)

  ggData[, gene := as.character(gene)]
  ggData[, AveExpr := as.numeric(AveExpr)]
  ggData[, logFC := as.numeric(logFC)]
  ggData[, adj.P.Val := as.numeric(adj.P.Val)]

  ggData[, neglog10FDR := -log10(adj.P.Val)]

  ggData[, gene_type := "ns"]
  ggData[!is.na(adj.P.Val) & !is.na(logFC) & adj.P.Val <= p_cut & logFC >= fc_up, gene_type := "up"]
  ggData[!is.na(adj.P.Val) & !is.na(logFC) & adj.P.Val <= p_cut & logFC <= fc_down, gene_type := "down"]

  ggData[, sig := !is.na(adj.P.Val) & adj.P.Val <= p_cut]

  # MA plot
  p_ma <- ggplot2::ggplot(ggData, ggplot2::aes(x = AveExpr, y = logFC, color = gene_type)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Pseudobulk MA",
      x = "Average expression",
      y = "logFC",
      color = "Class"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed")

  # Volcano plot
  p_vol <- ggplot2::ggplot(ggData, ggplot2::aes(x = logFC, y = neglog10FDR, color = gene_type)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Pseudobulk Volcano",
      x = "logFC",
      y = "-log10(FDR)",
      color = "Class"
    ) +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(fc_down, fc_up), linetype = "dashed")

  # Optional: label a few top hits (kept minimal and safe)
  top_lab <- ggData[order(adj.P.Val)][1:min(10, .N)]
  top_lab <- top_lab[!is.na(gene) & nzchar(gene) & is.finite(neglog10FDR) & is.finite(logFC)]
  if (nrow(top_lab) > 0 && sc_has_pkg("ggrepel")) {
    p_vol <- p_vol +
      ggrepel::geom_text_repel(
        data = top_lab,
        ggplot2::aes(label = gene),
        size = 5,                # increase font size
        fontface = "bold",       # optional but improves readability
        max.overlaps = 100,
        box.padding = 0.4,
        point.padding = 0.3,
        segment.size = 0.4
      )
  }

  list(
    p_ma = p_ma,
    p_vol = p_vol,
    tt = as.data.frame(ggData)
  )
}

############################################### UI ####################################################

scPseudobulk_ui <- function(id, sc1conf, sc1def) {
  
  ns <- NS(id)
  
  tabPanel(
    HTML("Pseudobulk DE (voom limma)"),
    h4("Pseudobulk differential expression using voom and limma"),
    "Reads raw counts from sc1counts.h5. Select replicate and condition, preview design and model, choose reference, then run DE.",
    br(), br(),
    
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        
        h4("Cell grouping"),
        selectInput(
          ns("sc1e1_celltype_col"),
          "Cell group column:",
          choices = sc1conf[grp == TRUE]$UI,
          selected = sc1conf[grp == TRUE]$UI[1]
        ),
        uiOutput(ns("sc1e1_celltype_level.ui")),
        
        br(),
        h4("Replicate and condition"),
        
        selectInput(
          ns("sc1e1_unit_col"),
          "Replicate column (sample or replicate):",
          choices = unique(c(sc1conf$UI)),
          selected = if ("orig.ident" %in% unique(c(sc1conf$UI))) "orig.ident" else unique(c(sc1conf$UI))[1]
        ),
        
        selectInput(
          ns("sc1e1_condition_col"),
          "Condition column (treatment):",
          choices = unique(c(sc1conf$UI)),
          selected = sc1conf[grp == TRUE]$UI[1]
        ),
        
        uiOutput(ns("sc1e1_reference_ui")),
        
        checkboxInput(
          ns("sc1e1_rep_within_cond"),
          "If replicate appears in multiple conditions, treat replicate within condition as a unique unit",
          value = FALSE
        ),
        
        br(),
        actionButton(ns("sc1e1_tog_cov"), "Optional covariates"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1e1_tog_cov")),
          selectInput(
            ns("sc1e1_covars"),
            "Covariates to include as columns in design data:",
            choices = unique(c(sc1conf$UI)),
            selected = character(0),
            multiple = TRUE
          )
        ),
        
        br(),
        actionButton(ns("sc1e1_tog_model"), "Model options"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1e1_tog_model")),
          textInput(
            ns("sc1e1_model"),
            "Design formula for model.matrix:",
            value = "~ condition"
          ),
          textInput(
            ns("sc1e1_contrast"),
            "Optional limma contrast (blank for default):",
            value = ""
          ),
          checkboxInput(ns("sc1e1_robust"), "Robust eBayes", value = TRUE)
        ),
        
        br(),
        actionButton(ns("sc1e1_tog_plot"), "Plot options"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1e1_tog_plot")),
          numericInput(ns("sc1e1_pcut"), "FDR cutoff", value = 0.05, step = 0.01),
          numericInput(ns("sc1e1_fc_up"), "Up logFC threshold", value = 1, step = 0.1),
          numericInput(ns("sc1e1_fc_down"), "Down logFC threshold", value = -1, step = 0.1),
          radioButtons(ns("sc1e1_fsz"), "Font size:", choices = c("Small","Medium","Large"),
                       selected = "Medium", inline = TRUE)
        ),
        
        br(),
        actionButton(ns("sc1e1_run"), "Run pseudobulk DE", class = "btn btn-primary")
      ),
      
      column(
        9,
        h4(uiOutput(ns("sc1e1_status"))),
        
        fluidRow(
          column(
            4,
            h4("Replicate by condition"),
            DT::DTOutput(ns("sc1e1_map"))
          ),
          column(
            4,
            h4("Design preview"),
            DT::DTOutput(ns("sc1e1_design_dt"))
          ),
          column(
            4,
            h4("Model matrix preview"),
            verbatimTextOutput(ns("sc1e1_model_txt"))
          )
        ),
        
        br(),
        
        h4("Pseudobulk Plots"),
        fluidRow(
          column(
            6,
            plotOutput(ns("sc1e1_ma"), height = "350px"),
            div(
              style = "display:flex; gap:8px; align-items:center;",
              downloadButton(ns("sc1e1_ma_pdf"), "MA PDF"),
              downloadButton(ns("sc1e1_ma_png"), "MA PNG")
            )
          ),
          column(
            6,
            plotOutput(ns("sc1e1_vol"), height = "450px"),
            div(
              style = "display:flex; gap:8px; align-items:center;",
              downloadButton(ns("sc1e1_vol_pdf"), "Volcano PDF"),
              downloadButton(ns("sc1e1_vol_png"), "Volcano PNG")
            )
          )
        ),
        
        br(),
        h4("Results table"),
        DT::DTOutput(ns("sc1e1_dt"))
      )
    )
  )
}

############################################### Server #################################################

scPseudobulk_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    observe_helpers()
    
    if (!exists("sList", inherits = TRUE)) {
      sList <<- c(Small = 10, Medium = 12, Large = 14)
    }
    
    h5_counts <- file.path(dir_inputs, "sc1counts.h5")
    
    output$sc1e1_celltype_level.ui <- renderUI({
      req(input$sc1e1_celltype_col)
      col <- sc_resolve_meta_col(sc1conf, sc1meta, input$sc1e1_celltype_col)
      req(col)
      
      lv <- as.character(sc1meta[[col]])
      lv <- lv[!is.na(lv) & nzchar(lv)]
      lv <- sort(unique(lv))
      req(length(lv) > 0)
      
      selectInput(
        ns("sc1e1_celltype_level"),
        "Cell group level:",
        choices = lv,
        selected = lv[1]
      )
    })
    
    filtered_cells_idx <- reactive({
      req(input$sc1e1_celltype_col, input$sc1e1_celltype_level)
      
      col <- sc_resolve_meta_col(sc1conf, sc1meta, input$sc1e1_celltype_col)
      req(col)
      
      ct <- sc1meta[[col]]
      which(!is.na(ct) & (as.character(ct) == as.character(input$sc1e1_celltype_level)))
    })
    
    output$sc1e1_reference_ui <- renderUI({
      keep <- filtered_cells_idx()
      req(length(keep) > 0)
      req(input$sc1e1_condition_col)
      
      cond_col <- sc_resolve_meta_col(sc1conf, sc1meta, input$sc1e1_condition_col)
      req(cond_col)
      
      cond_vec <- as.character(sc1meta[[cond_col]][keep])
      cond_levels <- sort(unique(cond_vec[!is.na(cond_vec) & nzchar(cond_vec)]))
      req(length(cond_levels) > 1)
      
      selectInput(
        ns("sc1e1_reference"),
        "Reference condition (control level):",
        choices = cond_levels,
        selected = cond_levels[1]
      )
    })
    
    preview_design <- reactive({
      keep <- filtered_cells_idx()
      req(length(keep) > 0)
      
      unit_col <- sc_resolve_meta_col(sc1conf, sc1meta, input$sc1e1_unit_col)
      cond_col <- sc_resolve_meta_col(sc1conf, sc1meta, input$sc1e1_condition_col)
      req(unit_col, cond_col)
      
      unit_vec <- as.character(sc1meta[[unit_col]][keep])
      cond_vec <- as.character(sc1meta[[cond_col]][keep])
      
      ok <- !is.na(unit_vec) & nzchar(unit_vec) & !is.na(cond_vec) & nzchar(cond_vec)
      unit_vec <- unit_vec[ok]
      cond_vec <- cond_vec[ok]
      
      unit_key <- unit_vec
      if (isTRUE(input$sc1e1_rep_within_cond)) {
        unit_key <- paste0(unit_vec, "___", cond_vec)
      }
      
      unit_condition <- tapply(cond_vec, unit_key, function(x) unique(x[!is.na(x) & nzchar(x)]))
      flag_multi <- vapply(unit_condition, function(x) length(x) != 1, logical(1))
      unit_condition_chr <- rep(NA_character_, length(unit_condition))
      unit_condition_chr[!flag_multi] <- vapply(unit_condition[!flag_multi], function(x) x[[1]], character(1))
      
      df <- data.frame(
        unit = names(unit_condition),
        condition = unit_condition_chr,
        stringsAsFactors = FALSE
      )
      
      covars <- character(0)
      if (!is.null(input$sc1e1_covars) && length(input$sc1e1_covars) > 0) covars <- input$sc1e1_covars
      if (length(covars) > 0) {
        covar_cols <- vapply(covars, function(z) sc_resolve_meta_col(sc1conf, sc1meta, z), character(1))
        covar_cols <- covar_cols[!is.na(covar_cols) & nzchar(covar_cols)]
        covar_cols <- unique(covar_cols)
        
        if (length(covar_cols) > 0) {
          covar_df <- lapply(covar_cols, function(col) {
            v <- as.character(sc1meta[[col]][keep][ok])
            tapply(v, unit_key, function(x) {
              ux <- unique(x[!is.na(x) & nzchar(x)])
              if (length(ux) == 0) NA_character_
              else if (length(ux) == 1) ux[[1]]
              else NA_character_
            })
          })
          covar_df <- as.data.frame(covar_df, stringsAsFactors = FALSE)
          names(covar_df) <- covar_cols
          covar_df$unit <- rownames(covar_df)
          rownames(covar_df) <- NULL
          
          df <- merge(df, covar_df, by = "unit", all.x = TRUE, sort = FALSE)
        }
      }
      
      df$multi_condition_flag <- flag_multi[match(df$unit, names(flag_multi))]
      
      df
    })
    
    output$sc1e1_map <- DT::renderDT({
      keep <- filtered_cells_idx()
      req(length(keep) > 0)
      
      unit_col <- sc_resolve_meta_col(sc1conf, sc1meta, input$sc1e1_unit_col)
      cond_col <- sc_resolve_meta_col(sc1conf, sc1meta, input$sc1e1_condition_col)
      req(unit_col, cond_col)
      
      unit_vec <- as.character(sc1meta[[unit_col]][keep])
      cond_vec <- as.character(sc1meta[[cond_col]][keep])
      
      tab <- as.data.frame.matrix(table(unit_vec, cond_vec, useNA = "ifany"))
      tab$replicate <- rownames(tab)
      tab <- tab[, c("replicate", setdiff(colnames(tab), "replicate")), drop = FALSE]
      
      DT::datatable(
        tab,
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE)
      )
    }, server = FALSE)
    
    output$sc1e1_design_dt <- DT::renderDT({
      df <- preview_design()
      req(nrow(df) > 0)
      
      DT::datatable(
        df,
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE)
      )
    }, server = FALSE)
    
    output$sc1e1_model_txt <- renderPrint({
      df <- preview_design()
      req(nrow(df) > 0)
      
      if (any(isTRUE(df$multi_condition_flag), na.rm = TRUE)) {
        cat("Some units map to multiple conditions.\n")
        cat("Fix replicate and condition columns, or enable replicate within condition.\n")
        return(invisible(NULL))
      }
      
      df2 <- df
      df2$condition <- factor(df2$condition)
      
      if (!is.null(input$sc1e1_reference) && nzchar(input$sc1e1_reference) && input$sc1e1_reference %in% levels(df2$condition)) {
        df2$condition <- stats::relevel(df2$condition, ref = input$sc1e1_reference)
      }
      
      fml <- stats::as.formula(input$sc1e1_model)
      mm <- stats::model.matrix(fml, data = df2)
      
      cat("Formula:\n")
      cat(input$sc1e1_model, "\n\n")
      cat("Reference:\n")
      cat(if (!is.null(input$sc1e1_reference)) input$sc1e1_reference else "None", "\n\n")
      cat("Model matrix columns:\n")
      cat(paste(colnames(mm), collapse = ", "), "\n\n")
      cat("Head of model matrix:\n")
      print(utils::head(mm, 10))
    })
    
    res <- eventReactive(input$sc1e1_run, {
      shiny::validate(shiny::need(file.exists(h5_counts), paste0(
        "Cannot find raw counts file: ", h5_counts, ". Create it first as sc1counts.h5 in dir_inputs."
      )))
      
      keep <- filtered_cells_idx()
      shiny::validate(shiny::need(length(keep) > 0, "No cells after grouping filter."))
      
      df_prev <- preview_design()
      shiny::validate(shiny::need(nrow(df_prev) > 1, "No valid pseudobulk units after preview."))
      
      if (any(isTRUE(df_prev$multi_condition_flag), na.rm = TRUE)) {
        shiny::validate(shiny::need(FALSE, "Some replicates map to multiple conditions. Check replicate and condition columns, or enable replicate within condition."))
      }
      
      covars <- character(0)
      if (!is.null(input$sc1e1_covars) && length(input$sc1e1_covars) > 0) covars <- input$sc1e1_covars
      
      pb <- sc_make_pseudobulk(
        sc1conf = sc1conf,
        sc1meta = sc1meta,
        sc1gene = sc1gene,
        h5_counts_path = h5_counts,
        keep_cells_idx = keep,
        unit_col_ui_or_meta = input$sc1e1_unit_col,
        condition_col_ui_or_meta = input$sc1e1_condition_col,
        covar_cols_ui_or_meta = covars,
        rep_within_cond = isTRUE(input$sc1e1_rep_within_cond),
        chunk_genes = 1000
      )
      
      fit <- sc_run_voom_limma(
        pseudobulk = pb$pseudobulk,
        design_df = pb$design_df,
        model_str = input$sc1e1_model,
        contrast_str = input$sc1e1_contrast,
        robust_ebayes = isTRUE(input$sc1e1_robust),
        reference_level = input$sc1e1_reference
      )
      
      plots <- sc_pseudobulk_plots(
        fit$table,
        p_cut = input$sc1e1_pcut,
        fc_up = input$sc1e1_fc_up,
        fc_down = input$sc1e1_fc_down
      )
      
      list(pb = pb, fit = fit, plots = plots)
    })
    
    output$sc1e1_status <- renderUI({
      if (!file.exists(h5_counts)) {
        return(HTML(paste0("Missing raw counts file: ", basename(h5_counts))))
      }
      if (is.null(res())) return(HTML("Ready"))
      HTML(paste0(
        "Done. Units: ", nrow(res()$pb$design_df),
        "   Genes: ", nrow(res()$pb$pseudobulk),
        "   Conditions: ", paste(levels(factor(res()$fit$design_df$condition)), collapse = ", "),
        "   Reference: ", input$sc1e1_reference
      ))
    })
    
    output$sc1e1_ma <- renderPlot({
      req(res())
      p <- res()$plots$p_ma +
        ggplot2::theme(text = ggplot2::element_text(size = sList[input$sc1e1_fsz]))
      print(p)
    })
    
    output$sc1e1_vol <- renderPlot({
      req(res())
      p <- res()$plots$p_vol +
        ggplot2::theme(text = ggplot2::element_text(size = sList[input$sc1e1_fsz]))
      print(p)
    })
    
    output$sc1e1_dt <- DT::renderDT({
      req(res())
      tt <- res()$plots$tt
      
      DT::datatable(
        tt,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          pageLength = 25,
          dom = "Bfrtip",
          buttons = c("excel", "csv"),
          scrollX = TRUE
        )
      )
    }, server = FALSE)
    
    output$sc1e1_ma_pdf <- downloadHandler(
      filename = function() {
        paste0("pseudobulk_MA_", input$sc1e1_celltype_level, ".pdf")
      },
      content = function(file) {
        req(res())
        p <- res()$plots$p_ma +
          ggplot2::theme(text = ggplot2::element_text(size = sList[input$sc1e1_fsz]))
        ggplot2::ggsave(file, plot = p, device = "pdf", width = 7, height = 5, useDingbats = FALSE)
      }
    )
    
    output$sc1e1_ma_png <- downloadHandler(
      filename = function() {
        paste0("pseudobulk_MA_", input$sc1e1_celltype_level, ".png")
      },
      content = function(file) {
        req(res())
        p <- res()$plots$p_ma +
          ggplot2::theme(text = ggplot2::element_text(size = sList[input$sc1e1_fsz]))
        ggplot2::ggsave(file, plot = p, device = "png", width = 7, height = 5, dpi = 300)
      }
    )
    
    output$sc1e1_vol_pdf <- downloadHandler(
      filename = function() {
        paste0("pseudobulk_Volcano_", input$sc1e1_celltype_level, ".pdf")
      },
      content = function(file) {
        req(res())
        p <- res()$plots$p_vol +
          ggplot2::theme(text = ggplot2::element_text(size = sList[input$sc1e1_fsz]))
        ggplot2::ggsave(file, plot = p, device = "pdf", width = 7, height = 6, useDingbats = FALSE)
      }
    )
    
    output$sc1e1_vol_png <- downloadHandler(
      filename = function() {
        paste0("pseudobulk_Volcano_", input$sc1e1_celltype_level, ".png")
      },
      content = function(file) {
        req(res())
        p <- res()$plots$p_vol +
          ggplot2::theme(text = ggplot2::element_text(size = sList[input$sc1e1_fsz]))
        ggplot2::ggsave(file, plot = p, device = "png", width = 7, height = 6, dpi = 300)
      }
    )
    
  })
}



############################################### Registration #################################################

register_tab(
  id     = "pseudobulk",
  title  = "Pseudobulk DE",
  ui     = scPseudobulk_ui,
  server = scPseudobulk_server
)