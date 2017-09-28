# R/3.4.1
library(shiny)
library(biomaRt)
library(ggplot2)

# Load essential data ----------------------------------------------------------
load("/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/RNA_seq_Rdata/RNA_seq_tpm.Rdata")
pheno <- read.delim("/projects/ytakemon/Col4a5xDO/Phenotype/Minimal_shiny_pheno.txt", sep = "\t", header = TRUE)
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "mmusculus_gene_ensembl",
                      verbose = TRUE)

# User Interface --------------------------------------------------------------
ui <- fluidPage(

  # Title
  tags$h2(tags$a("Col4a5 x Diversity Outbred", href = "http://ctshiny01:3838/ytakemon/Col4a5xDO/")," – Spearman Correlation: Phenotype v. Gene Expression"),

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
      # Select Sex
      selectInput(inputId = "sex",
                  label = "Select Sex",
                  choices = c("Both",
                              "Females",
                              "Males")),
      # Dowload eQTl map
      downloadButton("download_plot",
                     label = "Download"),
      br(),
      br(),
      div("Col4a5xDO Spearman Correlation: Phenotype v. Gene Expression v.1.0.0, powered by R/Shiny, developed by Yuka Takemon, ",
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
  Correlation_plot <- function(){
    gene_select <- input$gene_input # gene_select <- "Pik3r1"
    pheno_select <- input$pheno # pheno_select <- "ACR at 6 weeks"
    sex_select <- input$sex # sex_select <- "Males"

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

    # Update selected phenotype
    if (pheno_select == "Glomerular filtration rate"){
      pheno_select <- "C2_log"
    } else if (pheno_select == "ACR at 6 weeks"){
      pheno_select <- "ACR6WK_log"
    } else if (pheno_select == "ACR at 10 weeks"){
      pheno_select <- "ACR10WK_log"
    } else if (pheno_select == "ACR at 15 weeks"){
      pheno_select <- "ACR15WK_log"
    }

    #  Update selected Sex
    if (sex_select == "Both"){
      sex_select <- c("M","F")
    } else if (sex_select == "Females"){
      sex_select <- "F"
    } else if (sex_select == "Males"){
      sex_select <- "M"
    }

    # Extract selected gene TPM count
    TPM <- RNA_seq[,mart_extract$ensembl_gene_id]
    # Check to see if TPM was found
    if (nrow(TPM) == 0){
      validate(nrow(TPM) != 0, "Cannot query! Query gene not found in RNA-seq data.")
      stop("Cannot query! Query gene not found in RNA-seq data.")
    }
    pheno$tpm <- TPM

    # Subset pheno data based on selected sexes
    pheno_sub <- pheno[pheno$Sex %in% sex_select,]
    # Make sure complete cases for selected phenotype
    pheno_sub <- pheno_sub[complete.cases(pheno_sub[,pheno_select]),]

    # Calculate equation for subtitle
    fit <- lm(pheno_sub[,pheno_select] ~ pheno_sub[,"tpm"], data = pheno_sub)
    fitsum <- summary(fit)
    intcp <- signif(coef(fit)[1], 3)
    slope <- signif(coef(fit)[2], 3)
    pval <- signif(fitsum$coefficients[2,4], 3)
    r2 <- signif(fitsum$adj.r.squared, 3)
    eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")

# difference between cor and lm,r2
cor(pheno_sub[,pheno_select], pheno_sub[,"tpm"], method = "pearson")
###### TOBE CONTINUED









    ggplot(pheno_sub, aes(y = pheno_sub[, "tpm"], x = pheno_sub[, pheno_select])) +
          geom_smooth( method = lm) +
          geom_point() +
          scale_x_continuous(paste(gene_select, " TPM (Transcript per million)")) +
          scale_y_continuous(paste("Log-transformed ", input$pheno)) +
          labs( title = paste(""gene_select, " "))

data <- RNA_pheno[ complete.cases(RNA_pheno$C2_log),]

fit <- lm(C2_log ~ Tgm2, data)
fitsum <- summary(fit)
intcp <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(fitsum$coefficients[2,4], 3)
r2 <- signif(fitsum$adj.r.squared, 3)
eq <- paste("y = ", slope,"x ","+ ", intcp, ", ", "R^2 =", r2, ", ", " pval = ", pval, sep = "")
GFR <- ggplot(data, aes(y = C2_log, x = Tgm2)) +
	geom_smooth(method = lm) +
	geom_point() +
	background_grid(major = 'y', minor = "none") +
	annotate( "text" , x = 200, y = 6.5,label = eq, fontface = "bold", size = 4) +
	scale_x_continuous( "Tgm2 TPM", breaks = seq(0, max(data$Tgm2), by = 50)) +
	scale_y_continuous( "log-transformed C2 GFR", breaks = seq(0, max(data$C2_log), by = 0.5)) +
	labs( title = "Tgm2 TPM vs log C2 GFR (Both sexes, n=126)") +
	theme(plot.title = element_text(hjust = 0.5))










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
    Correlation_plot()
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
