library(shiny)
# All eQTL maps are located in helix @ /home/ytakemon/ShinyApps/Col4a5xDO/eQTL/www/

# User Interface --------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Col4a5 x Diversity Outbred – eQTL maps"),

  # Sidebar layout with input and output definitions ------------------
  sidebarLayout(
    # Sidebar panel for inputs ------------------
    sidebarPanel(
        p(span("Please use safari to view eQTL maps!", style = "color:blue")),
        p("Both Ensembl ID and gene names can be queried."),
        br(),
        textInput(inputId = "gene_input", label = "Gene Query", value = "Enter query here", width = "100%"),
        fluidRow(column(verbatimTextOutput("value"), width = 12))
        #fluidRow(column(widtextOutput("path_input")))
    ),
    # Main panel for displaying outputs ------------------
    mainPanel(
    imageOutput("image")
    )
  )
)
# Server ----------------------------------------------------------------------
server <- function(input, output) {
  # Display query gene in below text entry box
  output$value <- renderText({input$gene_input})

  # Render image
  output$image <- renderImage({
    # Find image
    path <- "/home/ytakemon/ShinyApps/Col4a5xDO/eQTL/www/"
    file <- list.files(path = path, pattern = input$gene_input)
    filename <- paste0(path,file)
    list(src = filename)
    }, deleteFile = FALSE)
}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
