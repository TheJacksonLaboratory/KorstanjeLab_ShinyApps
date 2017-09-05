# User Interface --------------------------------------------------------------
ui <- fluidPage(
  titlePanel("X-linked Alport Diversity Outbred mouse model – eQTL maps"),

  # Sidebar layout with input and output definitions ------------------
  sidebarLayout(
    # Sidebar panel for inputs ------------------
    sidebarPanel(
      textInput("Text", label = h3("Gene Query"), value = "Enter gene name or Ensembl ID"),
      fluidRow(column(3, verbatimTextOutput("value")))
    ),
    # Main panel for displaying outputs ------------------
    mainPanel(
      imageOutput("eQTL_map")
    )
  )
)

# Server ----------------------------------------------------------------------
server <- function(input, output) {
  output$value <- renderPrint({input$Text})

  image_path <-

}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)


Lrrc74a.eQT









setwd("/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/plot/RNA_qtl/Complete_eQTL_plots")
x <- list.files("/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/plot/RNA_qtl/Complete_eQTL_plots")
x[grep(paste0(".", "Lrrc74a", "."), x)]
