# Install all packages
install.packages(c("devtools", "shiny", "ggplot2", "reshape2","readxl", "boot","Matrix","mgcv"))
library(c("devtools", "shiny", "ggplot2", "reshape2","readxl", "boot","Matrix","mgcv"))
# Install bioconductor packages
library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite(c("qtl","biomaRt","DOQTL"))
