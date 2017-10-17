# Install all packages
install.packages(c("devtools", "shiny", "ggplot2", "reshape2","readxl", "boot","Matrix","mgcv"))
library(c("devtools", "shiny", "ggplot2", "reshape2","readxl", "boot","Matrix","mgcv"))
# Install bioconductor packages
library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite(c("qtl","biomaRt","DOQTL"))



# Shiny run logs are currently located:
# ytakemon@ctronshiny01:/var/log/shiny-server
# The apps can be found here:
# ytakemon@ctronshiny01:/opt/KorstanjeLab
