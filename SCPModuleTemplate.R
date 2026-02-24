#' Scaffold a new ShinyCell style module script
#'
#' Creates a minimal self contained module file with:
#' Functions block, UI block, Server block, Registration block.
#'
#' It also:
#' 1 checks existing module scripts in module_dir for ID conflicts
#' 2 picks the next available sc1a or sc1b or sc1c prefix number so namespaces are unique
#' 3 prevents filename conflicts
#'
#' Assumptions based on your module style
#' register_tab(id=..., title=..., ui=..., server=...)
#' UI input ids use ns("sc1xN...")
#' Server includes observe_helpers() and local fallbacks for pList pList2 pList3
#'
#' @param module_dir folder containing your module .R scripts
#' @param module_name base name for the file and functions, eg "scMyNewTab" (file will be scMyNewTab.R)
#' @param tab_id register_tab id, eg "bubble_heatmap"
#' @param tab_title register_tab title, eg "Bubble Plot"
#' @param tab_header UI title shown in tabPanel, eg "Bubbleplot / Heatmap"
#' @param prefix_group one of "a","b","c","d","e" to choose sc1a or sc1b etc
#' @param overwrite allow overwriting an existing file with same name
#' @return invisible list with file path and chosen prefix
scScaffoldModule <- function(module_dir,
                             module_name,
                             tab_id,
                             tab_title,
                             tab_header = tab_title,
                             prefix_group = c("a","b","c","d","e"),
                             overwrite = FALSE) {
  
  prefix_group <- match.arg(prefix_group)
  if (!dir.exists(module_dir)) stop("module_dir does not exist: ", module_dir)
  
  # list .R files
  r_files <- list.files(module_dir, pattern = "\\.R$", full.names = TRUE)
  all_text <- character(0)
  if (length(r_files) > 0) {
    all_text <- vapply(
      r_files,
      function(f) paste(readLines(f, warn = FALSE), collapse = "\n"),
      FUN.VALUE = character(1)
    )
  }
  
  # 1 check register_tab id conflict
  existing_tab_ids <- character(0)
  if (length(all_text) > 0) {
    m <- gregexpr("register_tab\\s*\\(.*?\\bid\\s*=\\s*['\\\"]([^'\\\"]+)['\\\"]", all_text, perl = TRUE)
    hits <- regmatches(all_text, m)
    if (length(hits) > 0) {
      existing_tab_ids <- unique(unlist(lapply(hits, function(x) {
        if (length(x) == 0) return(character(0))
        sub(".*\\bid\\s*=\\s*['\\\"]([^'\\\"]+)['\\\"].*", "\\1", x)
      })))
      existing_tab_ids <- existing_tab_ids[nzchar(existing_tab_ids)]
    }
  }
  if (tab_id %in% existing_tab_ids) {
    stop("tab_id conflicts with an existing register_tab id: ", tab_id)
  }
  
  # 2 pick next available sc1 prefix number for the chosen group
  # Search for any occurrence of sc1<group><number>
  used_nums <- integer(0)
  if (length(all_text) > 0) {
    pattern <- paste0("\\bsc1", prefix_group, "(\\d+)\\b")
    for (txt in all_text) {
      mm <- gregexpr(pattern, txt, perl = TRUE)
      hits <- regmatches(txt, mm)[[1]]
      if (length(hits) > 0) {
        nums <- suppressWarnings(as.integer(sub(paste0("^sc1", prefix_group), "", hits)))
        used_nums <- c(used_nums, nums[!is.na(nums)])
      }
    }
  }
  used_nums <- unique(used_nums)
  next_num <- if (length(used_nums) == 0) 1L else (max(used_nums) + 1L)
  
  ns_prefix <- paste0("sc1", prefix_group, next_num)  # eg sc1d2
  
  # 3 prevent filename conflicts
  out_file <- file.path(module_dir, paste0(module_name, ".R"))
  if (file.exists(out_file) && !isTRUE(overwrite)) {
    stop("Output file already exists: ", out_file, "\nSet overwrite = TRUE if you want to overwrite it.")
  }
  
  # function names follow your pattern: <module_name>_ui, <module_name>_server
  ui_fun <- paste0(module_name, "_ui")
  srv_fun <- paste0(module_name, "_server")
  
  # helper to build names like sc1d2inp1, sc1d2togL etc
  id_inp1 <- paste0(ns_prefix, "inp1")
  id_inp2 <- paste0(ns_prefix, "inp2")
  id_togL <- paste0(ns_prefix, "togL")
  id_sub1 <- paste0(ns_prefix, "sub1")
  id_sub2 <- paste0(ns_prefix, "sub2")
  id_oup  <- paste0(ns_prefix, "oup")
  id_psz  <- paste0(ns_prefix, "psz")
  id_fsz  <- paste0(ns_prefix, "fsz")
  id_oup_h <- paste0(id_oup, ".h")
  id_oup_w <- paste0(id_oup, ".w")
  
  # minimal template that matches your module structure and makes it self contained
  template <- c(
    paste0("# id     = \"", tab_id, "\""),
    paste0("# title  = \"", tab_title, "\""),
    "",
    "############################################### Functions ############################################",
    "",
    "# TODO implement your plot or table functions here",
    paste0("sc", ns_prefix, "_main <- function(inpConf, inpMeta, ...) {"),
    "  # return a ggplot object or a grob",
    "  stop(\"Not implemented\")",
    "}",
    "",
    "############################################### UI ####################################################",
    "",
    paste0(ui_fun, " <- function(id, sc1conf, sc1def) {"),
    "  ns <- NS(id)",
    "",
    "  tabPanel(",
    paste0("    HTML(\"", tab_header, "\"),"),
    paste0("    h4(\"", tab_header, "\"),"),
    "    br(), br(),",
    "",
    "    fluidRow(",
    "      column(",
    "        3, style = \"border-right: 2px solid black\",",
    paste0("        selectInput(ns(\"", id_inp1, "\"), \"Input 1:\", choices = sc1conf$UI, selected = sc1conf$UI[1]),"),
    paste0("        selectInput(ns(\"", id_inp2, "\"), \"Input 2:\", choices = sc1conf$UI, selected = sc1conf$UI[1]),"),
    "        br(),",
    paste0("        actionButton(ns(\"", id_togL, "\"), \"Filter Cells\"),"),
    "        conditionalPanel(",
    paste0("          condition = sprintf(\"input['%s'] %% 2 == 1\", ns(\"", id_togL, "\")),"),
    paste0("          selectInput(ns(\"", id_sub1, "\"), \"Cell information to subset:\", choices = sc1conf[grp == TRUE]$UI, selected = sc1def$grp1),"),
    paste0("          uiOutput(ns(\"", id_sub1, ".ui\")),"),
    paste0("          actionButton(ns(\"", id_sub1, "all\"), \"Select all groups\", class = \"btn btn-primary\"),"),
    paste0("          actionButton(ns(\"", id_sub1, "non\"), \"Deselect all groups\", class = \"btn btn-primary\")"),
    "        ),",
    "        br(), br(),",
    paste0("        actionButton(ns(\"", ns_prefix, "tog\"), \"Customize Plot\"),"),
    "        conditionalPanel(",
    paste0("          condition = sprintf(\"input['%s'] %% 2 == 1\", ns(\"", ns_prefix, "tog\")),"),
    paste0("          radioButtons(ns(\"", id_psz, "\"), \"Plot size:\", choices = c(\"Small\",\"Medium\",\"Large\"), selected = \"Medium\", inline = TRUE),"),
    paste0("          radioButtons(ns(\"", id_fsz, "\"), \"Font size:\", choices = c(\"Small\",\"Medium\",\"Large\"), selected = \"Medium\", inline = TRUE)"),
    "        )",
    "      ),",
    "",
    "      column(",
    paste0("        9, uiOutput(ns(\"", id_oup, ".ui\")),"),
    paste0("        downloadButton(ns(\"", id_oup, ".pdf\"), \"Download PDF\"),"),
    paste0("        downloadButton(ns(\"", id_oup, ".png\"), \"Download PNG\"),"),
    "        br(),",
    "        div(style = \"display:inline-block\",",
    paste0("            numericInput(ns(\"", id_oup_h, "\"), \"PDF / PNG height:\", width = \"138px\", min = 4, max = 20, value = 10, step = 0.5)"),
    "        ),",
    "        div(style = \"display:inline-block\",",
    paste0("            numericInput(ns(\"", id_oup_w, "\"), \"PDF / PNG width:\", width = \"138px\", min = 4, max = 20, value = 10, step = 0.5)"),
    "        )",
    "      )",
    "    )",
    "  )",
    "}",
    "",
    "############################################### Server #################################################",
    "",
    paste0(srv_fun, " <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {"),
    "  moduleServer(id, function(input, output, session) {",
    "    ns <- session$ns",
    "",
    "    # ShinyCell helper tags",
    "    observe_helpers()",
    "",
    "    # Self contained fallbacks for plot sizes",
    "    if (!exists(\"pList2\", inherits = TRUE)) {",
    "      pList2 <<- c(Small = \"350px\", Medium = \"550px\", Large = \"750px\")",
    "    }",
    "",
    "    # Filter UI",
    paste0("    output$", id_sub1, ".ui <- renderUI({"),
    paste0("      req(input$", id_sub1, ")"),
    paste0("      sub <- strsplit(sc1conf[UI == input$", id_sub1, "]$fID, \"\\\\|\")[[1]]"),
    "      checkboxGroupInput(",
    paste0("        ns(\"", id_sub2, "\"), \"Select which cells to show\","),
    "        inline = TRUE, choices = sub, selected = sub",
    "      )",
    "    })",
    "",
    paste0("    observeEvent(input$", id_sub1, "non, {"),
    paste0("      req(input$", id_sub1, ")"),
    paste0("      sub <- strsplit(sc1conf[UI == input$", id_sub1, "]$fID, \"\\\\|\")[[1]]"),
    "      updateCheckboxGroupInput(",
    "        session,",
    paste0("        inputId = \"", id_sub2, "\","),
    "        label = \"Select which cells to show\",",
    "        choices = sub, selected = NULL, inline = TRUE",
    "      )",
    "    })",
    "",
    paste0("    observeEvent(input$", id_sub1, "all, {"),
    paste0("      req(input$", id_sub1, ")"),
    paste0("      sub <- strsplit(sc1conf[UI == input$", id_sub1, "]$fID, \"\\\\|\")[[1]]"),
    "      updateCheckboxGroupInput(",
    "        session,",
    paste0("        inputId = \"", id_sub2, "\","),
    "        label = \"Select which cells to show\",",
    "        choices = sub, selected = sub, inline = TRUE",
    "      )",
    "    })",
    "",
    "    # Plot",
    paste0("    output$", id_oup, " <- renderPlot({"),
    paste0("      req(input$", id_inp1, ", input$", id_inp2, ")"),
    paste0("      sc", ns_prefix, "_main("),
    "        sc1conf, sc1meta,",
    paste0("        inp1 = input$", id_inp1, ","),
    paste0("        inp2 = input$", id_inp2, ","),
    paste0("        inpsub1 = input$", id_sub1, ","),
    paste0("        inpsub2 = input$", id_sub2, ","),
    paste0("        inpfsz  = input$", id_fsz, ","),
    "        dir_inputs = dir_inputs",
    "      )",
    "    })",
    "",
    paste0("    output$", id_oup, ".ui <- renderUI({"),
    paste0("      req(input$", id_psz, ")"),
    paste0("      plotOutput(ns(\"", id_oup, "\"), height = pList2[input$", id_psz, "])"),
    "    })",
    "",
    "    # Downloads",
    paste0("    output$", id_oup, ".pdf <- downloadHandler("),
    "      filename = function() {",
    paste0("        paste0(\"", ns_prefix, "_\", input$", id_inp1, ", \"_\", input$", id_inp2, ", \".pdf\")"),
    "      },",
    "      content = function(file) {",
    "        ggsave(",
    "          file, device = \"pdf\",",
    paste0("          height = input$", id_oup_h, ", width = input$", id_oup_w, ","),
    "          useDingbats = FALSE,",
    paste0("          plot = sc", ns_prefix, "_main("),
    "            sc1conf, sc1meta,",
    paste0("            inp1 = input$", id_inp1, ","),
    paste0("            inp2 = input$", id_inp2, ","),
    paste0("            inpsub1 = input$", id_sub1, ","),
    paste0("            inpsub2 = input$", id_sub2, ","),
    paste0("            inpfsz  = input$", id_fsz, ","),
    "            dir_inputs = dir_inputs",
    "          )",
    "        )",
    "      }",
    "    )",
    "",
    paste0("    output$", id_oup, ".png <- downloadHandler("),
    "      filename = function() {",
    paste0("        paste0(\"", ns_prefix, "_\", input$", id_inp1, ", \"_\", input$", id_inp2, ", \".png\")"),
    "      },",
    "      content = function(file) {",
    "        ggsave(",
    "          file, device = \"png\",",
    paste0("          height = input$", id_oup_h, ", width = input$", id_oup_w, ","),
    paste0("          plot = sc", ns_prefix, "_main("),
    "            sc1conf, sc1meta,",
    paste0("            inp1 = input$", id_inp1, ","),
    paste0("            inp2 = input$", id_inp2, ","),
    paste0("            inpsub1 = input$", id_sub1, ","),
    paste0("            inpsub2 = input$", id_sub2, ","),
    paste0("            inpfsz  = input$", id_fsz, ","),
    "            dir_inputs = dir_inputs",
    "          )",
    "        )",
    "      }",
    "    )",
    "",
    "  })",
    "}",
    "",
    "############################################### Registration #################################################",
    "",
    "register_tab(",
    paste0("  id     = \"", tab_id, "\","),
    paste0("  title  = \"", tab_title, "\","),
    paste0("  ui     = ", ui_fun, ","),
    paste0("  server = ", srv_fun),
    ")",
    ""
  )
  
  writeLines(template, out_file)
  
  invisible(list(
    file = out_file,
    prefix = ns_prefix,
    next_number = next_num,
    group = prefix_group
  ))
}