createSCPModuleTemplate <- function(
    module_dir,
    module_name,
    tab_id,
    tab_title,
    tab_header = tab_title,
    prefix_group = c("a","b","c","d","e"),
    overwrite = FALSE
) {
  message("createSCPModuleTemplate called")
  flush.console()
  
  prefix_group <- match.arg(prefix_group)
  if (!dir.exists(module_dir)) stop("module_dir does not exist: ", module_dir)
  
  r_files <- list.files(module_dir, pattern = "\\.R$", full.names = TRUE)
  
  all_text <- character(0)
  if (length(r_files) > 0) {
    all_text <- vapply(
      r_files,
      function(f) paste(readLines(f, warn = FALSE), collapse = "\n"),
      FUN.VALUE = character(1)
    )
  }
  
  ##########################################################################
  # 1 Check register_tab id conflicts
  ##########################################################################
  
  existing_tab_ids <- character(0)
  if (length(all_text) > 0) {
    pattern_id <- "register_tab\\s*\\(.*?\\bid\\s*=\\s*['\\\"]([^'\\\"]+)['\\\"]"
    m <- gregexpr(pattern_id, all_text, perl = TRUE)
    hits <- regmatches(all_text, m)
    
    if (length(hits) > 0) {
      existing_tab_ids <- unique(unlist(lapply(hits, function(x) {
        if (!length(x)) return(character(0))
        sub(pattern_id, "\\1", x, perl = TRUE)
      })))
      existing_tab_ids <- existing_tab_ids[nzchar(existing_tab_ids)]
    }
  }
  
  if (tab_id %in% existing_tab_ids) {
    stop("tab_id conflicts with an existing register_tab id: ", tab_id)
  }
  
  ##########################################################################
  # 2 Pick next available sc1<group><number>
  ##########################################################################
  
  used_nums <- integer(0)
  
  if (length(all_text) > 0) {
    # Match digits but ensure the next character is NOT another digit
    # This correctly matches sc1a2inp1 and sc1a2_main
    # and does NOT treat sc1a10 as sc1a1
    pattern_sc <- paste0("\\bsc1", prefix_group, "([0-9]+)(?![0-9])")
  

    for (txt in all_text) {
      m <- gregexpr(pattern_sc, txt, perl = TRUE)

      hits <- sort(unique(regmatches(txt, m)[[1]]))
      print("This are the matches from this modules")
      print(hits)
  
      
      if (length(hits) > 0) {
        nums <- suppressWarnings(as.integer(sub(pattern_sc, "\\1", hits, perl = TRUE)))
        used_nums <- c(used_nums, nums[!is.na(nums)])
       
      }
    }
  }
  
  used_nums <- unique(used_nums)
  next_num <- if (!length(used_nums)) 1L else (max(used_nums) + 1L)
  
  ns_prefix <- paste0("sc1", prefix_group, next_num)
  
  ##########################################################################
  # 3 Prevent filename conflicts
  ##########################################################################
  
  out_file <- file.path(module_dir, paste0(module_name, ".R"))
  if (file.exists(out_file) && !isTRUE(overwrite)) {
    stop("Output file already exists: ", out_file, "\nSet overwrite = TRUE if you want to overwrite it.")
  }
  
  ui_fun  <- paste0(module_name, "_ui")
  srv_fun <- paste0(module_name, "_server")
  
  id_inp1  <- paste0(ns_prefix, "inp1")
  id_inp2  <- paste0(ns_prefix, "inp2")
  id_togL  <- paste0(ns_prefix, "togL")
  id_sub1  <- paste0(ns_prefix, "sub1")
  id_sub2  <- paste0(ns_prefix, "sub2")
  id_oup   <- paste0(ns_prefix, "oup")
  id_psz   <- paste0(ns_prefix, "psz")
  id_fsz   <- paste0(ns_prefix, "fsz")
  id_oup_h <- paste0(id_oup, ".h")
  id_oup_w <- paste0(id_oup, ".w")
  
  template <- c(
    paste0("# id     = \"", tab_id, "\""),
    paste0("# title  = \"", tab_title, "\""),
    "",
    "############################################### Functions ############################################",
    "",
    "# TODO implement your plot or table functions here",
    # FIX: do not prepend an extra 'sc'
    paste0(ns_prefix, "_main <- function(inpConf, inpMeta, ...) {"),
    "  # return a ggplot object or a grob",
    "  stop(\"This is the template tab. Please add some code in the function section\")",
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
    "    observe_helpers()",
    "",
    "    if (!exists(\"pList2\", inherits = TRUE)) {",
    "      pList2 <<- c(Small = \"350px\", Medium = \"550px\", Large = \"750px\")",
    "    }",
    "",
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
    paste0("    output$", id_oup, " <- renderPlot({"),
    paste0("      req(input$", id_inp1, ", input$", id_inp2, ")"),
    # FIX: do not prepend extra 'sc'
    paste0("      ", ns_prefix, "_main("),
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
    paste0("    output$", id_oup, ".pdf <- downloadHandler("),
    "      filename = function() {",
    paste0("        paste0(\"", ns_prefix, "_\", input$", id_inp1, ", \"_\", input$", id_inp2, ", \".pdf\")"),
    "      },",
    "      content = function(file) {",
    "        ggsave(",
    "          file, device = \"pdf\",",
    paste0("          height = input$", id_oup_h, ", width = input$", id_oup_w, ","),
    "          useDingbats = FALSE,",
    paste0("          plot = ", ns_prefix, "_main("),
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
    paste0("          plot = ", ns_prefix, "_main("),
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
    group = prefix_group,
    used_numbers = sort(used_nums)
  ))
}