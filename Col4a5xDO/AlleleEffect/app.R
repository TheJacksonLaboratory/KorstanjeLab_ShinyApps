library(shiny)
library(biomaRt)
library(ggplot2)
library(DOQTL)
library(reshape2)

# Load essential data ---------------------------------------------------------
load("/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
load("/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/GM_snps.Rdata")
pheno <- read.delim("/projects/ytakemon/Col4a5xDO/Phenotype/Minimal_shiny_pheno.txt", sep = "\t", header = TRUE)


# Get mm10 data
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "mmusculus_gene_ensembl",
                      verbose = TRUE)
# User Interface --------------------------------------------------------------
ui <- fluidPage(

  # Title
  titlePanel("Col4a5 x Diversity Outbred – Allele Effect"),

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
      div("Col4a5xDO Allele Effect v.1.0.0, powered by R/Shiny, developed by Yuka Takemon, ",
          "souce code on ", a("Github", href = "https://github.com/TheJacksonLaboratory/KorstanjeLab_ShinyApps"),
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
  isoform_plot <- function(){
    gene_select <- input$gene_input
    pheno_select <- input$pheno

    #Select phenotype and load qtl file
    qtl_dir <- "/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/qtl/"
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
      stop("Cannot query! No marker found in query gene region.")
    }

    #calculate coef at target marker depending on the pheno_select variable
    fit <- lm(pheno$ACR6WK_log ~ pheno$Sex + best.genoprobs.192[,,target$marker], na.action = na.exclude)









































    if (is.na(names(gg_data)[2])){
      gg_data$Transcript <- transcript_list
      ggplot(gg_data, aes(x =Transcript, y = TPM, fill =Transcript)) +
      	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = binwidth, dotsize = dotsize) +
      	scale_x_discrete(paste0(gene, " Transcripts")) +
      	scale_y_continuous("TPM Counts") +
        guides(fill = FALSE)+
      	labs( title = paste0("Comparison of ", gene, " transcript TPM counts")) +
      	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 0, vjust = 0.5))
    } else {
    names(gg_data)[2] <- "Transcripts"
    ggplot(gg_data, aes(x =Transcripts, y = TPM, fill =Transcripts)) +
    	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = binwidth, dotsize = dotsize) +
    	scale_x_discrete(paste0(gene, " Transcripts")) +
    	scale_y_continuous("TPM Counts") +
      guides(fill = FALSE)+
    	labs( title = paste0("Comparison of ", gene, " transcript TPM counts")) +
    	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 15, vjust = 0.5))
    }
  }

  # Render plot ---------------------------------------
  output$plot <- renderPlot({
    isoform_plot()
    })

  # Download plot --------------------------
  output$download_plot <- downloadHandler(
    filename <- function(){
      paste0(input$gene_input, "_expressed_isoforms.pdf")
      },
    # content must be a function with arguemtn files to write plot
    content <- function(file) {
      pdf(file, width = 11, height = 7) #open device
        print(isoform_plot()) #print plot
      dev.off() # close device
    }
  )
}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
