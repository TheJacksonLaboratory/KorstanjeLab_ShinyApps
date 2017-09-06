library(shiny)
# User Interface --------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Col4a5 x Diversity Outbred – eQTL maps"),

  # Sidebar layout with input and output definitions ------------------
  sidebarLayout(
    # Sidebar panel for inputs ------------------
    sidebarPanel(
        textInput(inputId = "gene_input", label = "Gene Query", value = "Enter gene here", width = "100%"),
        fluidRow(column(verbatimTextOutput("value"), width = 12))
        #fluidRow(column(widtextOutput("path_input")))
    ),
    # Main panel for displaying outputs ------------------
    mainPanel( "main panel",
    img(src = "ENSMUSG00000020774.Aspa.eQTL.perm1000.pdf")
    )
  )
)
# Server ----------------------------------------------------------------------
server <- function(input, output) {
  # Display query gene in below text entry box
  output$value <- renderText({input$gene_input})

  # Render image
  path <- "/home/ytakemon/ShinyApps/Col4a5xDO/eQTL/www/"
  file <- list.files(path = path, pattern = "Aspa")
  output$image <- renderText({
    file
    })
  #output$image <- renderImage({
    #Find image
    #path <- "/home/ytakemon/ShinyApps/Col4a5xDO/eQTL/www/"
    #file <- list.files(path = path, pattern = "Aspa")
    #filename <- paste0(path,file)
    #list(src = "minion.jpeg")
    #}, deleteFile = FALSE)
}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
