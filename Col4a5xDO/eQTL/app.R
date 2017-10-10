library(shiny)
# All eQTL maps are located in helix @ /home/ytakemon/ShinyApps/Col4a5xDO/eQTL/www/
setwd("/home/ytakemon/ShinyApps/Col4a5xDO/eQTL/www/")
# User Interface --------------------------------------------------------------
ui <- fluidPage(

  # Title
  tags$h2(tags$a("Col4a5 x Diversity Outbred", href = "http://ctshiny01:3838/ytakemon/Col4a5xDO/")," – eQTL"),

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
      div("Col4a5xDO eQTL v.1.0.1, powered by R/Shiny, developed by ",
          a("Yuka Takemon", href="mailto:yuka.takemon@jax.org?subject=KorstanejeLab shiny page"),
          ", souce code on ", a("Github", href = "https://github.com/TheJacksonLaboratory/KorstanjeLab_ShinyApps"),
          " (JAX network only).")
    ),

    # Main panel for displaying outputs ------------------
    mainPanel(
      imageOutput("image")
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
  file <- list.files(pattern = paste0("\\.",input$gene_input,"\\."))

  if(input$gene_input == "Enter query here"){
    validate(
      need(input$gene_input != "Enter query here", "Please enter a gene.")
    )
  } else if (input$gene_input != "Enter query here"){
    validate(
      need(length(file) == 1, "Queried gene cannot be found.")
    )
  }


  list(src = file)
  }, deleteFile = FALSE)

  # Download handler --------------------------
  output$download_image <- downloadHandler(
    filename <- paste0(input$gene_input, "_eQTL_map.pdf"),
    content <- function(downloadFile) {
      file <- list.files(pattern = input$gene_input)
      file.copy(file, downloadFile)
    }
  )
}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
