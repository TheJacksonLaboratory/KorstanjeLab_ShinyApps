library(shiny)
library(biomaRt)
library(ggplot2)
library(reshape2)

# Load essential data ---------------------------------------------------------
load("/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/Gene_allele.Rdata")
Gene_allele <- as.data.frame(Gene_allele)
load("/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/RNA_seq_Rdata/All_transcript_tpm.Rdata")
# Get mm10 data
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "mmusculus_gene_ensembl",
                      verbose = TRUE)
# User Interface --------------------------------------------------------------
ui <- fluidPage(

  # Title
  titlePanel("Col4a5 x Diversity Outbred – Isoform query"),

  # Sidebar layout with input and output definitions ------------------
  sidebarLayout(
    # Sidebar panel for inputs ------------------
    sidebarPanel(
      # Header
      #p(span("Please use safari to view eQTL maps!", style = "color:blue")),
      p("Both Ensembl ID and gene names can be queried."),
      br(),
      # Input text box
      textInput(inputId = "gene_input",
                label = "Gene Query",
                value = "Enter query here",
                width = "100%"),
      # Confirmation text
      fluidRow(column(verbatimTextOutput("value"),
               width = 12)),
      # Slider Input for bin width
      sliderInput(inputId = "binwidth",
                  label = "Adjust bin width:",
                  min = 0,
                  max = 1,
                  value = 0.02,
                  step = 0.02),
      # Slider input for dot size
      sliderInput(inputId = "dotsize",
                  label = "Adjust dot size:",
                  min = 0,
                  max = 10,
                  value = 0.2,
                  step = 0.2),
      # Dowload eQTl map
      downloadButton("download_plot",
                     label = "Download"),
      br(),
      br(),
      div("Col4a5xDO Isoforms v.1.0.0, powered by R/Shiny, developed by Yuka Takemon, ",
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
    gene <- input$gene_input
    binwidth <- input$binwidth
    dotsize <- input$dotsize

    # Figure out if gene symbol or ENSEMBL ID
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
    # Extract transcript list for input gene
    transcript_list <- mart_extract$ensembl_transcript_id
    gene_id <- mart_extract$ensembl_gene_id[1]

    #Create new dataframe with only TPM from list of transcirpts
    transcript_list <- colnames(All_transcript_tpm)[colnames(All_transcript_tpm) %in% transcript_list]
    TPM <- All_transcript_tpm[,transcript_list]

    #plot non-transformed TPM
    gg_data <- melt(TPM, id.vars = transcript_list,
    	value.name = "TPM")

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
