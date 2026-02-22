############################################## Functions ################################################
  

if(!is.na((marker_list_roc_path)) && file.exists(marker_list_roc_path) && !dir.exists(marker_list_roc_path)){
top_genes_selection <- c('10','20','30','50','all')

############################################## UI ################################################


scMarkerGenes_ui <- function(id, sc1conf, sc1def){
    tabPanel( 
      HTML("Markers Lists ROC"), 
      h4("markers lists ROC"), 
      "In this tab, users can subset the markers list per cluster given each resolution calculated using ROC ",  
      br(),br(), 
      fluidRow(
        sidebarLayout(
          # sidebarPanel(
          inputPanel(
            selectInput("resolution","Clustering resolution:",
                        choices= sc1conf$UI[grep("res",sc1conf$UI)],
                        selected  = sc1conf$UI[grep("res",sc1conf$UI)][1],multiple = FALSE),
            selectInput("assay","Select assay:",choices= assays, selected  = assays[1]),
            selectInput("top","Number of genes per cluster:",choices= top_genes_selection, selected  = top_genes_selection[1])
            
          ),
          mainPanel(
            DTOutput("mytable_markers_roc")
            
          )
        )  
        
      )
    )
}
 
############################################## Server ################################################

 # this can be modified to use the output of presto::wilcoxauc directly. 
scBubbHeat_server <- function(id, sc1conf, sc1def) {
    ### ????? added table ROC
    output$mytable_markers_roc <- renderDT(server = FALSE,{
      resolution_assay<- paste0(input$resolution)
      top_selection <- input$top
      
      if(top_selection != 'all'){ 
        
        top_gene<-as.numeric(top_selection)
        DT::datatable(
          markers_list_roc %>% filter(res  == resolution_assay) %>% group_by(cluster) %>% top_n(n = top_gene , wt = log2FC), 
          extensions = c("Buttons"),
          options = list(
            dom = 'Bfrtip',
            buttons = list(
              list(extend = "csv", text = "Download Top Genes", filename =paste0("SubsetMarkersList_",resolution_assay),
                  exportOptions = list(
                    modifier = list(page = "all")
                  )
              )
            )
          )
        )
      }else{
        DT::datatable(
          markers_list_roc %>% filter(res  == resolution_assay) %>% group_by(cluster), 
          extensions = c("Buttons"),
          options = list(
            dom = 'Bfrtip',
            buttons = list(
              list(extend = "csv", text = "Download Full Results", filename =paste0("MarkersList_",resolution_assay),
                  exportOptions = list(
                    modifier = list(page = "all")
                  )
              )
            )
          )
        )
      }
      
    })
  } 
 
  