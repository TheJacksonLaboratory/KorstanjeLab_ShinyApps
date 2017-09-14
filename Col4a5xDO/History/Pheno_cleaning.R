# Clean up pheno data to make Shiny work faster

# read and clean up phenotype data
load("/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/genoprobs/best.genoprobs.192.Rdata")
pheno <- read.delim("/projects/ytakemon/Col4a5xDO/Phenotype/1415_master_pheno.txt",
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE)

rownames(pheno) <- make.names(pheno[,1]) #move sample ID to row names
pheno <- pheno[rownames(best.genoprobs.192),] #subset pheno to match 192
#clean up pheno and add log of ACR
pheno[pheno < 0 ] = NA
pheno[pheno ==  -Inf] = NA
pheno$ACR6WK_log <- log(pheno$ACR6WK)
pheno$ACR10WK_log <- log(pheno$ACR10WK)
pheno$ACR15WK_log <- log(pheno$ACR15WK)
pheno$GFRSigma <- as.numeric(gsub(",","",pheno[,"GFRSigma"]))
pheno$C2[81] = NA #over sigma 650000 cut off
pheno$GFRSigma[81] = NA #over sigma 650000 cut off
pheno$C2[148] = NA #over sigma 650000 cut off
pheno$GFRSigma[148] = NA #over sigma 650000 cut off
pheno$C2_log <- log(pheno$C2)
pheno[pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs
# Subset only needed data
pheno <- pheno[,c("MouseID", "Sex", "C2", "C2_log", "ACR6WK", "ACR6WK_log",
                  "ACR10WK", "ACR10WK_log", "ACR15WK", "ACR15WK_log")]

write.table(pheno, "/projects/ytakemon/Col4a5xDO/Phenotype/Minimal_shiny_pheno.txt",
            sep = "\t")
