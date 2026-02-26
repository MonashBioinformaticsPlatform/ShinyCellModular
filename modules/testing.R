# id     = "testing"
# title  = "testing"

############################################### Functions ############################################

# TODO implement your plot or table functions here
sc1a4_main <- function(inpConf, inpMeta, ...) {
  # return a ggplot object or a grob
  stop("Not implemented")
}

############################################### UI ####################################################

testing_ui <- function(id, sc1conf, sc1def) {
  ns <- NS(id)

  tabPanel(
    HTML("testing"),
    h4("testing"),
    br(), br(),

    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        selectInput(ns("sc1a4inp1"), "Input 1:", choices = sc1conf$UI, selected = sc1conf$UI[1]),
        selectInput(ns("sc1a4inp2"), "Input 2:", choices = sc1conf$UI, selected = sc1conf$UI[1]),
        br(),
        actionButton(ns("sc1a4togL"), "Filter Cells"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a4togL")),
          selectInput(ns("sc1a4sub1"), "Cell information to subset:", choices = sc1conf[grp == TRUE]$UI, selected = sc1def$grp1),
          uiOutput(ns("sc1a4sub1.ui")),
          actionButton(ns("sc1a4sub1all"), "Select all groups", class = "btn btn-primary"),
          actionButton(ns("sc1a4sub1non"), "Deselect all groups", class = "btn btn-primary")
        ),
        br(), br(),
        actionButton(ns("sc1a4tog"), "Customize Plot"),
        conditionalPanel(
          condition = sprintf("input['%s'] %% 2 == 1", ns("sc1a4tog")),
          radioButtons(ns("sc1a4psz"), "Plot size:", choices = c("Small","Medium","Large"), selected = "Medium", inline = TRUE),
          radioButtons(ns("sc1a4fsz"), "Font size:", choices = c("Small","Medium","Large"), selected = "Medium", inline = TRUE)
        )
      ),

      column(
        9, uiOutput(ns("sc1a4oup.ui")),
        downloadButton(ns("sc1a4oup.pdf"), "Download PDF"),
        downloadButton(ns("sc1a4oup.png"), "Download PNG"),
        br(),
        div(style = "display:inline-block",
            numericInput(ns("sc1a4oup.h"), "PDF / PNG height:", width = "138px", min = 4, max = 20, value = 10, step = 0.5)
        ),
        div(style = "display:inline-block",
            numericInput(ns("sc1a4oup.w"), "PDF / PNG width:", width = "138px", min = 4, max = 20, value = 10, step = 0.5)
        )
      )
    )
  )
}

############################################### Server #################################################

testing_server <- function(id, sc1conf, sc1meta, sc1gene, sc1def, dir_inputs) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observe_helpers()

    if (!exists("pList2", inherits = TRUE)) {
      pList2 <<- c(Small = "350px", Medium = "550px", Large = "750px")
    }

    output$sc1a4sub1.ui <- renderUI({
      req(input$sc1a4sub1)
      sub <- strsplit(sc1conf[UI == input$sc1a4sub1]$fID, "\\|")[[1]]
      checkboxGroupInput(
        ns("sc1a4sub2"), "Select which cells to show",
        inline = TRUE, choices = sub, selected = sub
      )
    })

    observeEvent(input$sc1a4sub1non, {
      req(input$sc1a4sub1)
      sub <- strsplit(sc1conf[UI == input$sc1a4sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(
        session,
        inputId = "sc1a4sub2",
        label = "Select which cells to show",
        choices = sub, selected = NULL, inline = TRUE
      )
    })

    observeEvent(input$sc1a4sub1all, {
      req(input$sc1a4sub1)
      sub <- strsplit(sc1conf[UI == input$sc1a4sub1]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(
        session,
        inputId = "sc1a4sub2",
        label = "Select which cells to show",
        choices = sub, selected = sub, inline = TRUE
      )
    })

    output$sc1a4oup <- renderPlot({
      req(input$sc1a4inp1, input$sc1a4inp2)
      sc1a4_main(
        sc1conf, sc1meta,
        inp1 = input$sc1a4inp1,
        inp2 = input$sc1a4inp2,
        inpsub1 = input$sc1a4sub1,
        inpsub2 = input$sc1a4sub2,
        inpfsz  = input$sc1a4fsz,
        dir_inputs = dir_inputs
      )
    })

    output$sc1a4oup.ui <- renderUI({
      req(input$sc1a4psz)
      plotOutput(ns("sc1a4oup"), height = pList2[input$sc1a4psz])
    })

    output$sc1a4oup.pdf <- downloadHandler(
      filename = function() {
        paste0("sc1a4_", input$sc1a4inp1, "_", input$sc1a4inp2, ".pdf")
      },
      content = function(file) {
        ggsave(
          file, device = "pdf",
          height = input$sc1a4oup.h, width = input$sc1a4oup.w,
          useDingbats = FALSE,
          plot = sc1a4_main(
            sc1conf, sc1meta,
            inp1 = input$sc1a4inp1,
            inp2 = input$sc1a4inp2,
            inpsub1 = input$sc1a4sub1,
            inpsub2 = input$sc1a4sub2,
            inpfsz  = input$sc1a4fsz,
            dir_inputs = dir_inputs
          )
        )
      }
    )

    output$sc1a4oup.png <- downloadHandler(
      filename = function() {
        paste0("sc1a4_", input$sc1a4inp1, "_", input$sc1a4inp2, ".png")
      },
      content = function(file) {
        ggsave(
          file, device = "png",
          height = input$sc1a4oup.h, width = input$sc1a4oup.w,
          plot = sc1a4_main(
            sc1conf, sc1meta,
            inp1 = input$sc1a4inp1,
            inp2 = input$sc1a4inp2,
            inpsub1 = input$sc1a4sub1,
            inpsub2 = input$sc1a4sub2,
            inpfsz  = input$sc1a4fsz,
            dir_inputs = dir_inputs
          )
        )
      }
    )

  })
}

############################################### Registration #################################################

register_tab(
  id     = "testing",
  title  = "testing",
  ui     = testing_ui,
  server = testing_server
)

