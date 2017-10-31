library(ggplot2)
library(biomaRt)
library(shiny)
load("/opt/KorstanjeLab/Internal/AgingDO/data/shiny_annotation.RData")
load("/opt/KorstanjeLab/Internal/AgingDO/data/DO188b_kidney.RData")
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "mmusculus_gene_ensembl",
                      verbose = TRUE)
# User Interface --------------------------------------------------------------
ui <- fluidPage(

  # Title
  tags$h2(tags$a("Aging DO", href = "./Internal/AgingDO/")," – Age * QTL Interaction"),

  # Sidebar layout with input and output definitions ------------------
  sidebarLayout(
    # Sidebar panel for inputs ------------------
    sidebarPanel(
      # Header
      p("Only MGI symbols can be queried."),
      br(),
      # Select level
      selectInput(inputId = "level",
                  label = "Select level",
                  choices = c("mRNA",
                              "protein")),
      # Input text box
      textInput(inputId = "gene_input",
                label = "Gene Query",
                value = "Enter query here",
                width = "100%"),
      # Confirmation text
      fluidRow(column(verbatimTextOutput("value"),
               width = 12)),
      # Select chromosome
      selectInput(inputId = "chr",
                  label = "Select chromosome",
                  choices = c(1:19,"X")),
      # Dowload eQTl map
      downloadButton("download_plot",
                     label = "Download"),
      br(),
      br(),
      div("Aging DO: Age * QTL Interaction v.1.0.0, powered by R/Shiny, maintained by ",
          a("Yuka Takemon", href="mailto:yuka.takemon@jax.org?subject=KorstanejeLab shiny page"),
          ", souce code on ",
          a("Github", href = "https://github.com/TheJacksonLaboratory/KorstanjeLab_ShinyApps", target = "_blank"),
          " (JAX network only).",
          br(),
          "Connect with us @",
          a("The Korstanje Lab", href = "https://www.jax.org/research-and-faculty/research-labs/the-korstanje-lab", target = "_blank"))
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
  AgeQTLIntPlot <- function(){
    gene_select <- input$gene_input # gene_select <- "Lcat"
    level_select <- input$level # level_select <- "protein"
    chr_select <- input$chr # chr_select <- 15

    if (gene_select == "Enter query here"){
      return(NULL)
    }

    # for given MGI symbol, find Ensembl ids
    other.ids <- function(gene.name, level) {
      if (level == "mRNA") {
        sel <- which(mRNA.list$symbol == gene.name)[1]
        if (!is.na(sel)) return(mRNA.list[sel,]) else return(c(NA,NA,NA))
      }
      if (level == "protein") {
        sel <- which(protein.list$symbol == gene.name)[1]
        if (!is.na(sel)) return(protein.list[sel,]) else return(c(NA,NA,NA))
      }
    }

    # Validate
    query <- other.ids(gene_select, level_select)
    validate(
      need(nrow(query) == 1, "Query is not found in RNAseq or Proteomics dataset.")
    )

    # Gather gene info for query
    if (level_select == "mRNA"){
      id <- query$id[1]
      symbol <- query$symbol[1]
      file <- paste0("/opt/KorstanjeLab/Internal/AgingDO/data/Intscan_mrna/Age/",
                     id,
                     "_",
                     symbol,
                     ".rds")
    } else {
      id <- query$protein_id[1]
      symbol <- query$symbol[1]
      file <- paste0("/opt/KorstanjeLab/Internal/AgingDO/data/Intscan_prot/Age/",
                     id,
                     "_",
                     symbol,
                     ".rds")
    }

    # Validate file existance
    if (file.exists(file)){
      fit <- readRDS(file)
    } else {
      fit <- NA
    }
    validate(
      need(!is.na(fit), "Query is not found in RNAseq or Proteomics dataset.")
    )

    # fit
    fit <- as.data.frame(fit)
    fit$chr <- substr(rownames(fit),1,regexpr("_",rownames(fit)) - 1)
    fit.chr <- fit[fit$chr==chr_select,]
    max.marker <- rownames(fit.chr)[which.max(fit.chr$pheno1)]
    if (level_select=="mRNA") ens <- other.ids(symbol, level_select)[[1]] else ens <- other.ids(symbol, level_select)[[3]]
    # Validate existance before moving forward
    validate(
      need(!is.na(ens), "Query is not found in RNAseq or Proteomics dataset.")
    )
    if (level_select=="mRNA") y <- expr.mrna[,ens] else y <- expr.protein[,ens]
    y[annot.samples$Sex=="M"] <- y[annot.samples$Sex=="M"] - mean(y[annot.samples$Sex=="M"], na.rm=TRUE)
    y[annot.samples$Sex=="F"] <- y[annot.samples$Sex=="F"] - mean(y[annot.samples$Sex=="F"], na.rm=TRUE)
    coef <- se <- tstat <-  NULL

    # Get effect by age
    for (age in c("6", "12", "18")) {
      sel <- annot.samples$Age == age
      lm.fit <- summary(lm(y[sel] ~ 0 + genoprobs[sel,,max.marker]))
      coef <- cbind(coef, lm.fit$coef[1:8,1])
      se   <- cbind(se, lm.fit$coef[1:8,2])
      tstat   <- cbind(tstat, lm.fit$coef[1:8,3])
    }

    # Compile into df
    dt <- data.frame(Allele = LETTERS[rep(1:8,3)],
               Age = factor(rep(c(6,12,18), each=8)),
               beta = as.vector(coef),
               beta.se = as.vector(se),
               Tstat=as.vector(tstat))

    ggplot(aes(x=as.numeric(Allele)+as.numeric(Age)/15-2/15, y=beta,colour=Age), data=dt) +
          geom_errorbar(aes(ymin=beta-beta.se, ymax=beta+beta.se), width=.1) +
          geom_point(aes(size=abs(Tstat))) + xlab("Allele") +
          ylab("beta +/- SE") +
          scale_x_continuous(breaks=c(1:8), labels=LETTERS[1:8]) +
          geom_abline(intercept = 0, slope = 0, colour=I("grey")) +
          ggtitle(paste(max.marker, gene_select)) +
          theme(panel.border = element_rect(linetype = "solid", fill=NA, colour = "grey"))
  }

  # Render plot ---------------------------------------
  output$plot <- renderPlot({
      if (is.null(AgeQTLIntPlot())){
        return(NULL)
      }
      AgeQTLIntPlot()
    })

  # Download plot --------------------------
  output$download_plot <- downloadHandler(
    filename <- function(){
      paste0(input$gene_input, "_chr", input$chr,"_ageQTLInteraction.pdf")
      },
    # content must be a function with arguemtn files to write plot
    content <- function(file) {
      pdf(file, width = 11, height = 7) #open device
        print(AgeQTLIntPlot()) #print plot
      dev.off() # close device
    }
  )

  # Output links --------------------------------
  output$gene_links <- renderText({

    # if not ready
    if(is.null(AgeQTLIntPlot())){
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
