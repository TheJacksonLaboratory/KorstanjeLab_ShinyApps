library(shiny)
# All eQTL maps are located in helix @ /home/ytakemon/ShinyApps/Col4a5xDO/phenoQTL/www/
setwd("/home/ytakemon/ShinyApps/Col4a5xDO/phenoQTL/www/")

# User Interface --------------------------------------------------------------
ui <- fluidPage(

  # Title
  tags$h2(tags$a("Col4a5 x Diversity Outbred", href = "http://ctshiny01:3838/ytakemon/Col4a5xDO/")," – Phenotype QTL"),

  # Sidebar layout with input and output definitions ------------------
  sidebarLayout(
    # Sidebar panel for inputs ------------------
    sidebarPanel(
      # Header
      p(span("Please use safari to view phenoQTL maps!", style = "color:blue")),
      br(),
      # Dropdown menu
      selectInput(inputId = "pheno",
                  label = "Select Phenotype",
                  choices = c("Glomerular filtration rate",
                              "Albumin normalized to creatinine at 6 weeks",
                              "Albumin normalized to creatinine at 10 weeks",
                              "Albumin normalized to creatinine at 15 weeks")),
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

  # Render image -----------------------------
  output$image <- renderImage({
    # Find image
    # Based on input drop down choices
    if (input$pheno == "Glomerular filtration rate"){
      file <- "Figure4.1_qtl.log.C2.GFR.noX.pdf"
    } else if (input$pheno == "Albumin normalized to creatinine at 6 weeks"){
      file <- "Figure4.2_qtl.log.Alb6wk.noX.pdf"
    } else if (input$pheno == "Albumin normalized to creatinine at 10 weeks"){
      file <- "Figure4.3_qtl.log.Alb10wk.noX.pdf"
    } else if (input$pheno == "Albumin normalized to creatinine at 15 weeks"){
      file <- "Figure4.4_qtl.log.Alb15wk.noX.pdf"
    }
    # Define absolute path and list source
    list(src = file)
  }, deleteFile = FALSE)

  # Download handler --------------------------
  output$download_image <- downloadHandler(
    filename <- paste0(input$pheno, "_QTL_map.pdf"),
    content <- function(downloadFile) {
      # Based on input drop down choices
      if (input$pheno == "Glomerular filtration rate"){
        file <- "Figure4.1_qtl.log.C2.GFR.noX.pdf"
      } else if (input$pheno == "Albumin normalized to creatinine at 6 weeks"){
        file <- "Figure4.2_qtl.log.Alb6wk.noX.pdf"
      } else if (input$pheno == "Albumin normalized to creatinine at 10 weeks"){
        file <- "Figure4.3_qtl.log.Alb10wk.noX.pdf"
      } else if (input$pheno == "Albumin normalized to creatinine at 15 weeks"){
        file <- "Figure4.4_qtl.log.Alb15wk.noX.pdf"
      }
      file.copy(file, downloadFile)
    }
  )
}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
