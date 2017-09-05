# This is a trial to make the 1415 Col4a5xDO project data more accessible.
# Will be maintained by Yuka Takemon @ytakemon

# REQUIREMENTS:
# This shinyapp will be hosted via helix
# Users will need access to both helix and ctshiny01 servers
# Using R/3.4.1
# RStudio Server: http://ctshiny01:8787
# Shiny App Server: http://ctshiny01:3838
# Images must be in a directory call "www"

# NOTES:
# Any application that you put in your ~/ShinyApps directory is automatically hosted.
# The directory must be on helix or ctshiny01.
# The URL will be in the format  http://ctshiny01:3838/[YOUR USERNAME]/[APP NAME]
# http://ctshiny01:3838/ytakemon/Shiny_Col4a5xDO

# RESOURCES:
# layout guide: https://shiny.rstudio.com/articles/layout-guide.html

# DATA TRANSFER FROM CADILLAC TO HELIX
ssh cadillac
cd /hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct
scp -r best.compiled.genoprob ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO
# prompt will apepar to enter data


# LIBRARY INSTALLATION
# install.packages("devtools")
# source("http://bioconductor.org/biocLite.R")
# biocLite("DOQTL")
# install.packages("ggplot2")
# install.packages("reshape2")

setwd("/home/ytakemon/ShinyApps/")
library(ggplot2)
library(DOQTL)
library(reshape2)
library(shiny)
library(rmakrdown)

runApp("Shiny_Col4a5xDO")
