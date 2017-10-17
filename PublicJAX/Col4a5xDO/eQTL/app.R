library(shiny)
library(biomaRt)

# All eQTL maps are located in helix @ /home/ytakemon/ShinyApps/Col4a5xDO/eQTL/www/
setwd("/opt/KorstanjeLab/Col4a5xDO/eQTL/www/")
# User Interface --------------------------------------------------------------
ui <- fluidPage(

  # Title
  tags$h2(tags$a("Col4a5 x Diversity Outbred", href = "http://ctronshiny01:3838/KorstanjeLab/Col4a5xDO")," – eQTL"),

  # Sidebar layout with input and output definitions ------------------
  sidebarLayout(
    # Sidebar panel for inputs ------------------
    sidebarPanel(
      # Header
      p(span("Please use safari to view eQTL maps!", style = "color:blue")),
      p("Both Ensembl ID and gene names can be queried."),
      br(),
      # Input text box
      textInput(inputId = "gene_input", label = "Gene Query", value = "Enter query here", width = "100%"),
      # Confirmation text
      fluidRow(column(verbatimTextOutput("value"), width = 12)),
      # Dowload eQTl map
      downloadButton("download_image", label = "Download"),

      br(),
      br(),
      div("Col4a5xDO eQTL v.1.2.0, powered by R/Shiny, developed and maintained by ",
          a("Yuka Takemon", href="mailto:yuka.takemon@jax.org?subject=KorstanejeLab shiny page"),
          ", souce code on ",
          a("Github", href = "https://github.com/TheJacksonLaboratory/KorstanjeLab_ShinyApps", target = "_blank"),
          " (JAX network only).",
          br(),
          "Connect with us @",
          a("The Korstanje Lab", href = "https://www.jax.org/research-and-faculty/research-labs/the-korstanje-lab", target = "_blank"))
    ),

    # Main panel for displaying outputs ------------------
    mainPanel(
      imageOutput("image"),
      htmlOutput("gene_links")
    )
  )
)
# Server ----------------------------------------------------------------------
server <- function(input, output) {

  # Display query gene in below text entry box -----------------------
  output$value <- renderText({input$gene_input})


  # Render image -----------------------------
  output$image <- renderImage({
  # Find image
  file <- list.files(pattern = paste0("\\.",input$gene_input,"\\.")) # file <- list.files(pattern = paste0("\\.","Aspa","\\."))
  # Validate
  if(input$gene_input == "Enter query here"){
    validate(
      need(input$gene_input != "Enter query here", "Please enter a gene.")
    )
  } else if (input$gene_input != "Enter query here"){
    validate(
      need(length(file) == 1, "Queried gene cannot be found.")
    )
  }

  # Find page dimensions:
  # Source file
  list( src = file)
  }, deleteFile = FALSE)

  # Download handler --------------------------
  output$download_image <- downloadHandler(
    filename <- paste0(input$gene_input, "_eQTL_map.pdf"),
    content <- function(downloadFile) {
      file <- list.files(pattern = input$gene_input)
      file.copy(file, downloadFile)
    }
  )

  # Output links --------------------------------
  output$gene_links <- renderText({
    # if not ready
    file <- list.files(pattern = paste0("\\.",input$gene_input,"\\."))
    # Validate
    if(input$gene_input == "Enter query here"){
      validate(
        need(input$gene_input != "Enter query here", "Please enter a gene.")
      )
    } else if (input$gene_input != "Enter query here"){
      validate(
        need(length(file) == 1, "Queried gene cannot be found.")
      )
    }
    if(!file.exists(file)){
      return(NULL)
    }

    # find gene links
    gene <- input$gene_input

    ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                          dataset = "mmusculus_gene_ensembl",
                          verbose = TRUE)

    if ( substr(gene, 1, 7) == "ENSMUSG"){
      mart_extract <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "ensembl_transcript_id",
                                    "chromosome_name", "start_position", "end_position"),
                                    filters = "ensembl_gene_id",
                                    values = gene,
                                    mart = ensembl)
                                  } else {
      mart_extract <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "ensembl_transcript_id",
                                    "chromosome_name", "start_position", "end_position"),
                                    filters = "mgi_symbol",
                                    values = gene,
                                    mart = ensembl)
      }
    symbol <- mart_extract$mgi_symbol[1]
    ens_id <- mart_extract$ensembl_gene_id[1]
    chr <- mart_extract$chromosome_name[1]
    start <- mart_extract$start_position[1]
    end <- mart_extract$end_position[1]
    ensembl_link <- paste0("http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=", ens_id)
    mgi_link <- paste0("http://www.informatics.jax.org/searchtool/Search.do?query=", ens_id)

  # output links
  paste(br(),
        br(),
        br(),
        br(),
        br(),
        br(),
        br(),
        p(symbol, "is located on chromosome", chr, ":", start, "-", end),
        a("[Ensembl]", href = ensembl_link, target="_blank"),
        a("[MGI]", href = mgi_link, target = "_blank"))
})

}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
