# Condensed server.r and ui.r for Jax website
library(shiny)
library(readxl)
library(ggplot2)

source("/opt/KorstanjeLab/GFRcalc/helpers.R")
options(stringsAsFactors = FALSE)
max_plots <- 16

# User Interface --------------------------------------------------------------
ui <- pageWithSidebar(

  # Application title
  headerPanel("Glomerular Filtration Rate Estimation"),

  # Sidebar with a slider input for number of observations
  sidebarPanel(
    fileInput("GFRfile", "Input file:"),
    helpText("Please follow the template provided in this example:"),
    downloadButton("download_tempalate",
                   label = "Download Example Input File"),
    br(),
    br(),
    sliderInput("dilution","Dilution:", min=50, max=200, value=100),
    checkboxInput("Animal", "Animal (A-Z)", value=TRUE),
    checkboxInput("MouseId", "Mouse ID", value=TRUE),
    checkboxInput("Weight", "Weight", value=TRUE),
    checkboxInput("InjectedVolume", "Injected volume", value=TRUE),
    checkboxInput("GFR2", "Two-compartment model (2C)", value=TRUE),
    checkboxInput("GFR4", "2 One-compartment models (2xC1)", value=TRUE),
    checkboxInput("GFR1", "One-compartment model (1C)", value=TRUE),
    checkboxInput("GFR3", "Piecewise linear model (PL)", value=TRUE),
    checkboxInput("Sigma.C2", "Square root of the error variance (Sigma.C2)", value=TRUE),
    checkboxInput("nNA", "Number of missing observations (nNA)", value=TRUE),
  #  sliderInput("trhold","Outlier threshold:", min=1, max=10, value=5),
    br(),
    div("GFRcalc 0.9.5, powered by R/Shiny, developed by ",
        a("Petr Simecek", href = "https://github.com/simecek", target = "_blank"),
        " and currently maintained by ",
        a("Yuka Takemon", href="mailto:yuka.takemon@jax.org?subject=KorstanejeLab shiny page"), ". ",
        "Source code on ",
        a("Github", href="https://github.com/simecek/GFRcalc", target = "_blank"),".",
        br(),
        "Connect with us @",
        a("The Korstanje Lab", href = "https://www.jax.org/research-and-faculty/research-labs/the-korstanje-lab", target = "_blank")),
        div(a(img(src = "JAXlogo_trans.gif", height= 98, width= 186), href="https://www.jax.org"))
  ),

  mainPanel(
    textOutput("check"),
    tableOutput('contents'),
    uiOutput("plots")
  )
 )
# Server ----------------------------------------------------------------------
server <- function(input, output){

  # read data file
  file.upload <- reactive({
    inFile <- input$GFRfile
    if (!is.null(inFile)) { # if any file name is given
      ext <- tools::file_ext(inFile$name)
      file.rename(inFile$datapath, paste(inFile$datapath, ext, sep=".")) # dirty hack, see issue 85

      output <- read_excel(paste(inFile$datapath, ext, sep="."), sheet=2, col_names = FALSE)
      # output <- read_excel(paste(inFile, sep = "."), sheet=2, col_names = FALSE)
      output <- output[!is.na(output[,1]),] # remove NA rows
      output <- as.data.frame(output)
      if (!is.character(output[,1])) { # if animal IDs missing, add them
        animalID <- rep(LETTERS, each=length(unique(output[,1])))[1:nrow(output)]
        output <- cbind(animalID, output)
      }
      output <- output[,1:5] # only 5 columns matters
      names(output) <- c("Animal", "Time", "M1", "M2", "M3")

      # adding information about animals - should be on sheet 3
      animal.table <- read_excel(paste(inFile$datapath, ext, sep="."), sheet=3, col_names = FALSE)[2:4,-1]
      # animal.table <- read_excel(paste(inFile, sep = "."), sheet=3, col_names = FALSE)[2:4,-1]
      animal.table <- t(animal.table)
      attr(output, "animals") <- animal.table

      output
    } else {
      NULL
    }
  })

  output$check <- renderText({dt <- file.upload(); check.format(dt)})

  output$contents <- renderTable({

    dt <- file.upload()

    if (!is.null(dt) & check.format(dt)=="") { # everything ok

      # info about animals given? if not add them

      animals <- unique(dt$Animal)
      animal.table <- attr(dt, "animals")
      results <- NULL

      for (a in animals) {

        tmp <- subset(dt, Animal==a)
        tmp <- tmp[order(tmp$Time),]

        # very rough outlier detection
        tmp[-1,c("M1","M2","M3")][outlier.detect(tmp[-1,], 5)] <- NA

        tmp$mean <- rowMeans(tmp[,c("M1","M2","M3")], na.rm=TRUE)
        start <- 2 # skip first observation
        tmp2 <- tmp[start:nrow(tmp),]

        fit1 <- oneexp(y=tmp2$mean,x=tmp2$Time)
        fit2 <- twoexp(y=tmp2$mean,x=tmp2$Time)
        fit3 <- linint(y=tmp2$mean,x=tmp2$Time)
        fit4 <- twoexp.approx(y=tmp2$mean,x=tmp2$Time)

        inj.volume <- as.numeric(animal.table[animals==a,3])

        newrow <- data.frame(Animal = a,
                             MouseId = animal.table[animals==a,1],
                             Weight = as.numeric(animal.table[animals==a,2]),
                             "Injected Volume" = inj.volume,

                             "C2" = input$dilution*inj.volume*tmp$mean[1]*ifelse(is.null(fit2), NA, coef(fit2)[2]*coef(fit2)[4]/(coef(fit2)[1]*coef(fit2)[4]+coef(fit2)[2]*coef(fit2)[3])),
                             "2xC1" = ifelse(is.null(fit4),NA,input$dilution*inj.volume*tmp$mean[1]/fit4),
                             "C1" = input$dilution*inj.volume*tmp$mean[1]*ifelse(is.null(fit1), NA, coef(fit1)[2]/coef(fit1)[1]),
                             "PL" = ifelse(is.null(fit3), NA, input$dilution*inj.volume*tmp$mean[1]/fit3),
                             "Sigma.C2" = ifelse(is.null(fit2),NA,as.integer(round(summary(fit2)$sigma))),
                             nNA = sum(is.na(tmp2[,c("M1","M2","M3")])),
                             check.names = FALSE)

        results <- rbind(results, newrow)

      }

      if (!input$Animal) results$"Animal" <- NULL
      if (!input$MouseId) results$"MouseId" <- NULL
      if (!input$Weight) results$"Weight" <- NULL
      if (!input$InjectedVolume) results$"Injected Volume" <- NULL
      if (!input$GFR2) results$"C2" <- NULL
      if (!input$GFR4) results$"2xC1" <- NULL
      if (!input$GFR1) results$"C1" <- NULL
      if (!input$GFR3) results$"PL" <- NULL
      if (!input$Sigma.C2) results$"Sigma.C2" <- NULL
      if (!input$nNA) results$"nNA" <- NULL

      rownames(results) <- NULL

      return(results)

    } else {
      return (NULL)
    }
  }, digits=1, include.rownames=FALSE)

  output$plots <- renderUI({

    dt <- file.upload()

    if (!is.null(dt) & check.format(dt)=="") { # everything ok

      plot_output_list <- lapply(1:length(unique(dt[,1])), function(i) {
        plotname <- paste0("plot", i)
        plotOutput(plotname, height = 390, width = 520)
      })

      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, plot_output_list)
    }
  })

  for (i in 1:max_plots){
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste0("plot", my_i)
      output[[plotname]] <- renderPlot({make.plot(file.upload(), my_i)})
    })
  }

  # Download template ---------------------------------------
  # https://stackoverflow.com/questions/39449544/shiny-download-an-excel-file
  output$download_tempalate <- downloadHandler(
    filename = function(){
      paste("GFR_sample_file.xlsx")
    },
    content = function(file){
      myfile <-  "/opt/KorstanjeLab/GFRcalc/GFR_template.xlsx"
      file.copy(myfile, file)
    }
  )

}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
