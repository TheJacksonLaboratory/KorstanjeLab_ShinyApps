library(shiny)
library(biomaRt)
library(ggplot2)
library(reshape2)

# Load essential data ---------------------------------------------------------
setwd("/opt/KorstanjeLab/Col4a5xDO/RefData/")
load("./Gene_allele.Rdata")
Gene_allele <- as.data.frame(Gene_allele)
load("./All_transcript_tpm.Rdata")
# Get mm10 data
ensembl <- readRDS("./ensembl.rds")

# User Interface --------------------------------------------------------------
ui <- fluidPage(

  # Title
  tags$h2(tags$a("Col4a5 x Diversity Outbred", href = "/KorstanjeLab/Col4a5xDO/")," – Isoform query"),

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
                  max = 10,
                  value = 1,
                  step = 0.02),
      # Slider input for dot size
      sliderInput(inputId = "dotsize",
                  label = "Adjust dot size:",
                  min = 0,
                  max = 10,
                  value = 1,
                  step = 0.2),
      # Checkbox for allele plot
      checkboxInput(inputId = "allele_check",
                    label = "Plot by allele",
                    value = FALSE),

      # Dowload eQTl map
      downloadButton("download_plot",
                     label = "Download"),
      br(),
      br(),
      div("Col4a5xDO Isoforms v.1.2.0, powered by R/Shiny, developed and maintained by ",
          a("Yuka Takemon", href="mailto:yuka.takemon@jax.org?subject=KorstanejeLab shiny page"),
          ", souce code on ",
          a("Github", href = "https://github.com/TheJacksonLaboratory/KorstanjeLab_ShinyApps", target = "_blank"),
          " (JAX network only).",
          br(),
          "Connect with us @",
          a("The Korstanje Lab", href = "https://www.jax.org/research-and-faculty/research-labs/the-korstanje-lab", target = "_blank")),
          div(a(img(src = "JAXlogo_trans.gif", height= 98, width= 186), href="https://www.jax.org"))
    ),

    # Main panel for displaying outputs ------------------
    mainPanel(
      imageOutput("plot"),
      htmlOutput("gene_links")
    )
  )
)

# Server ----------------------------------------------------------------------
server <- function(input, output) {

  # Display query gene in below text entry box -----------------------
  output$value <- renderText({input$gene_input})

  # Create plot function ------------------------------------------
  isoform_plot <- function(){
    gene <- input$gene_input # gene <- "Aspa"
    binwidth <- input$binwidth # binwidth <- 0.02
    dotsize <- input$dotsize # dotsize <- 0.02

    if (gene == "Enter query here"){
      return(NULL)
    }


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

    # Validate query
    validate(
      need(nrow(mart_extract) > 0, "Gene not found in Ensembl mm10 database. Please check and try again!")
    )

    # Extract transcript list for input gene
    transcript_list <- mart_extract$ensembl_transcript_id
    gene_id <- mart_extract$ensembl_gene_id[1]

    #Create new dataframe with only TPM from list of transcirpts
    transcript_list <- colnames(All_transcript_tpm)[colnames(All_transcript_tpm) %in% transcript_list]
    TPM <- All_transcript_tpm[,transcript_list]
    gg_data <- melt(TPM, id.vars = transcript_list,
    	value.name = "TPM")
    gg_data$Allele <- Gene_allele[,mart_extract$ensembl_gene_id[1]]

    # Determine if allele check (default is FALSE)
    # If FALSE, allele check is not selected
    if (!input$allele_check){
    # iso_plot
      if (is.na(names(gg_data)[2])){
        gg_data$Transcript <- transcript_list
        ggplot(gg_data, aes(x =Transcript, y = TPM, fill =Transcript)) +
        	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = binwidth, dotsize = dotsize) +
          scale_fill_brewer(palette = "Set3") +
        	scale_x_discrete(paste0(gene, " Transcripts")) +
        	scale_y_continuous("TPM Counts") +
          guides(fill = FALSE)+
        	labs( title = paste0("Comparison of ", gene, " transcript TPM counts")) +
        	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 0, vjust = 0.5))
      } else {
      names(gg_data)[2] <- "Transcripts"
      ggplot(gg_data, aes(x =Transcripts, y = TPM, fill =Transcripts)) +
      	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = binwidth, dotsize = dotsize) +
        scale_fill_brewer(palette = "Set3") +
      	scale_x_discrete(paste0(gene, " Transcripts")) +
      	scale_y_continuous("TPM Counts") +
        guides(fill = FALSE)+
      	labs( title = paste0("Comparison of ", gene, " transcript TPM counts")) +
      	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 15, vjust = 0.5))
      }
    # If TRUE, allele check is on
    } else if (input$allele_check){
    # allele_iso_plot
      if (is.na(names(gg_data)[2])){
        gg_data$Transcript <- transcript_list
        ggplot(gg_data, aes(x =Transcript, y = TPM, fill =Allele)) +
        	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = binwidth, dotsize = dotsize, position_dodge(1)) +
        	scale_x_discrete(paste0(gene, " Transcripts")) +
        	scale_y_continuous("TPM Counts") +
          scale_fill_brewer(palette = "Set3") +
        	labs( title = paste0("Comparison of ", gene, " transcript TPM counts by allele")) +
        	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 0, vjust = 0.5))
      } else {
      names(gg_data)[2] <- "Transcripts"
      ggplot(gg_data, aes(x =Transcripts, y = TPM, fill =Allele)) +
      	geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2, dotsize = dotsize, position = position_dodge(1)) +
      	scale_x_discrete(paste0(gene, " Transcripts")) +
      	scale_y_continuous("TPM Counts") +
        scale_fill_brewer(palette = "Set3") +
      	labs( title = paste0("Comparison of ", gene, " transcript TPM counts by allele")) +
      	theme( plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 15, vjust = 0.5))
      }
    }
  }

  # Render plot ---------------------------------------
  output$plot <- renderPlot({
      if (is.null(isoform_plot())){
        return(NULL)
      }
      isoform_plot()
    })

  # Download plot -----------------------------------
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

  # Output links ----------------------------------
  output$gene_links <- renderText({
    # if not ready
    if(is.null(isoform_plot())){
      return(NULL)
    }

    # find gene links
    gene <- input$gene_input
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
    symbol <- mart_extract$mgi_symbol[1]
    ens_id <- mart_extract$ensembl_gene_id[1]
    chr <- mart_extract$chromosome_name[1]
    start <- mart_extract$start_position[1]
    end <- mart_extract$end_position[1]
    ensembl_link <- paste0("http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=", ens_id)
    mgi_link <- paste0("http://www.informatics.jax.org/searchtool/Search.do?query=", ens_id)

    # output links
    paste(p(symbol, "is located on chromosome", chr, ":", start, "-", end),
          a("[Ensembl]", href = ensembl_link, target="_blank"),
          a("[MGI]", href = mgi_link, target = "_blank"),br(),
          a("Diversity Outbred founder strain guide", href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4602074/figure/Fig1/", target = "_blank"))
  })

}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
