install.packages("devtools")
install.packages(c("devtools", "shiny", "ggplot2", "reshape2", "biomaRt","readxl"))

source("https://bioconductor.org/biocLite.R")
biocLite(c("DOQTL", "qtl"))
