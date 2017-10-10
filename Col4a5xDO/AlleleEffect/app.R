library(shiny)
library(biomaRt)
library(ggplot2)
library(DOQTL)
library(reshape2)

# Load essential data ---------------------------------------------------------
setwd("/home/ytakemon/ShinyApps/Col4a5xDO/RefData/")
load("best.genoprobs.192.Rdata")
load("GM_snps.Rdata")
pheno <- read.delim("Minimal_shiny_pheno.txt", sep = "\t", header = TRUE)

# Get mm10 data
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "mmusculus_gene_ensembl",
                      verbose = TRUE)
# User Interface --------------------------------------------------------------
ui <- fluidPage(

  # Title
  tags$h2(tags$a("Col4a5 x Diversity Outbred", href = "http://ctshiny01:3838/ytakemon/Col4a5xDO/")," – Allele Effect"),

  # Sidebar layout with input and output definitions ------------------
  sidebarLayout(
    # Sidebar panel for inputs ------------------
    sidebarPanel(
      # Header
      p("Both Ensembl ID and gene names can be queried."),
      br(),
      # Select Phenotype
      selectInput(inputId = "pheno",
                  label = "Select Phenotype",
                  choices = c("Glomerular filtration rate",
                              "ACR at 6 weeks",
                              "ACR at 10 weeks",
                              "ACR at 15 weeks")),
      # Input text box
      textInput(inputId = "gene_input",
                label = "Gene Query",
                value = "Enter query here",
                width = "100%"),
      # Confirmation text
      fluidRow(column(verbatimTextOutput("value"),
               width = 12)),
      # Dowload eQTl map
      downloadButton("download_plot",
                     label = "Download"),
      br(),
      br(),
      div("Col4a5xDO Allele Effect v.1.1.1, powered by R/Shiny, developed by ",
          a("Yuka Takemon", href="mailto:yuka.takemon@jax.org?subject=KorstanejeLab shiny page"),
          ", souce code on ", a("Github", href = "https://github.com/TheJacksonLaboratory/KorstanjeLab_ShinyApps"),
          " (JAX network only).")
    ),

    # Main panel for displaying outputs ------------------
    mainPanel(
      imageOutput("plot")
    )
  )
)

# Server ----------------------------------------------------------------------
server <- function(input, output) {

  # Display query gene in below text entry box -----------------------
  output$value <- renderText({input$gene_input})

  # Create plot function ------------------------------------------
  AlleleEffect_plot <- function(){
    gene_select <- input$gene_input # gene_select <- "Pik3r1"
    pheno_select <- input$pheno # pheno_select <- "ACR at 6 weeks"

    if (gene_select == "Enter query here"){
      return(NULL)
    }

    #Select phenotype and load qtl file
    qtl_dir <- "./qtl/"
    if (pheno_select == "Glomerular filtration rate"){
      file <- "qtl.GFR.log.C2.192.Rdata"
      load(paste0(qtl_dir, file))
      qtl <- qtl.GFR.log.C2.192
    } else if (pheno_select == "ACR at 6 weeks"){
      file <- "qtl.log.ACR6WK.192.Rdata"
      load(paste0(qtl_dir, file))
      qtl <- qtl.log.ACR6WK.192
    } else if (pheno_select == "ACR at 10 weeks"){
      file <- "qtl.log.ACR10WK.192.Rdata"
      load(paste0(qtl_dir, file))
      qtl  <- qtl.log.ACR10WK.192
    } else if (pheno_select == "ACR at 15 weeks"){
      file <- "qtl.log.ACR15WK.192.Rdata"
      load(paste0(qtl_dir, file))
      qtl <- qtl.log.ACR15WK.192
    }

    # Figure out if gene input is gene symbol or ENSEMBL_ID
    if ( substr(gene_select, 1, 7) == "ENSMUSG"){
      mart_extract <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name",
                                           "start_position", "end_position"),
                            filters = "ensembl_gene_id",
                            values = gene_select,
      											mart = ensembl)
    } else {
      mart_extract <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name",
                                           "start_position", "end_position"),
                            filters = "mgi_symbol",
                            values = gene_select,
      											mart = ensembl)
    }
    # Gather chr, start, end information
    chr <- mart_extract$chromosome_name
    start <- (mart_extract$start_position / (1e6))
    end <- (mart_extract$end_position / (1e6))

    # Find marker in QTL object
    qtl <- qtl$lod$A
    qtl <- qtl[qtl$chr == chr,]
    qtl <- qtl[qtl$pos >= start,]
    qtl <- qtl[qtl$pos <= end,]
    target <- qtl[qtl$lod == max(qtl$lod),]

    # Check to see if target was found
    if (nrow(target) == 0){
      validate(nrow(target) != 0, "Cannot query! No markers found in query gene region")
      stop("Cannot query! No marker found in query gene region.")
    }

    #calculate coef at target marker depending on the pheno_select variable
    if (pheno_select == "Glomerular filtration rate"){
      fit <- lm(pheno$C2_log ~ pheno$Sex + best.genoprobs.192[,,target$marker], na.action = na.exclude)
    } else if (pheno_select == "ACR at 6 weeks"){
      fit <- lm(pheno$ACR6WK_log ~ pheno$Sex + best.genoprobs.192[,,target$marker], na.action = na.exclude)
    } else if (pheno_select == "ACR at 10 weeks"){
      fit <- lm(pheno$ACR10WK_log ~ pheno$Sex + best.genoprobs.192[,,target$marker], na.action = na.exclude)
    } else if (pheno_select == "ACR at 15 weeks"){
      fit <- lm(pheno$ACR15WK_log ~ pheno$Sex + best.genoprobs.192[,,target$marker], na.action = na.exclude)
    }

    #Female averages by strain at Snp marker
    cfit <- coef(fit)
    cfit[10] = 0 #set standard (ie strain compared to) to 0
    cfit <- cfit + cfit[1]  #add female
    ecfit <- exp(cfit) #exp fixes values back to natural numbers
    ecfit_F <- as.data.frame(ecfit)

    #Male averages by strain at snp marker
    cfit <- coef(fit)
    cfit[10] = 0 #set standard (ie strain compared to) to 0
    cfit <- cfit + cfit[1] + cfit[2] #add female and male = male (female = female, male = male + female)
    ecfit <- exp(cfit) #exp fixes balues back to natural numbers
    ecfit_M <- as.data.frame(ecfit)

    # Combine both sexes
    data_coef <- ecfit_F
    data_coef$M <- ecfit_M$ecfit
    colnames(data_coef) <- c("F","M")
    names(data_coef) <- c("F","M")
    data_coef <- data_coef[c(3:10),]
    rownames(data_coef) <- c("A","B","C","D","E","F","G","H")
    data_coef$founder <- rownames(data_coef)

    ggdata <- melt(data_coef)
    colnames(ggdata) <- c("Founder", "Sex", "Value")
    ggplot(ggdata, aes( x = Founder, y = Value, fill = Sex)) +
    	geom_bar( stat = "identity", position = "dodge") +
    	scale_fill_discrete( labels = c("Females","Males")) +
    	labs(title = paste0("Col4a5xDO allele effect at ", gene_select ," (Chr", target$chr, " ", target$marker, " position: ", target$pos, ") for ",pheno_select, " by founder strains"),
    				x = "DO Founders",
    				y = paste0(pheno_select, " values")) +
    	theme( legend.position = "right", plot.title = element_text(hjust = 0.5))
    }

  # Render plot ---------------------------------------
  output$plot <- renderPlot({
      if (is.null(AlleleEffect_plot())){
        return(NULL)
      }
      AlleleEffect_plot()
    })

  # Download plot --------------------------
  output$download_plot <- downloadHandler(
    filename <- function(){
      paste0(input$gene_input, "_allele_effect.pdf")
      },
    # content must be a function with arguemtn files to write plot
    content <- function(file) {
      pdf(file, width = 11, height = 7) #open device
        print(AlleleEffect_plot()) #print plot
      dev.off() # close device
    }
  )
}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
