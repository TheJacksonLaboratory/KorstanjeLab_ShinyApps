# This is a trial to make the 1415 Col4a5xDO project data more accessible.
# Will be maintained by Yuka Takemon @ytakemon

# REQUIREMENTS:
# This shinyapp will be hosted via helix
# Users will need access to both helix and ctshiny01 servers
# Using R/3.4.1
# RStudio Server: http://ctshiny01:8787
# Shiny App Server: http://ctshiny01:3838

# NOTES:
# Any application that you put in your ~/ShinyApps directory is automatically hosted.
# The directory must be on helix or ctshiny01.
# The URL will be in the format  http://ctshiny01:3838/[YOUR USERNAME]/[APP NAME]

# LIBRARY INSTALLATION
# install.packages("devtools")
# source("http://bioconductor.org/biocLite.R")
# biocLite("DOQTL")
library(ggplot2)
library(DOQTL)
setwd("/projects/Shiny_Col4a5xDO")
