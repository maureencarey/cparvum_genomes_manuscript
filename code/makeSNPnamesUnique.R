# About -------------------------------------------------------------------
# 17/05/2018
#
# Reference: https://www.biostars.org/p/315219/
#
# 1000G plink bim file has duplicated IDs.
# - read plink BIM file
# - make unique SNPids
# - output BIM file
#
# input: BIM filname
# Note: This script will overwrite the existing input bim file


# Input arguments ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

fileBim <- args[1]

# Data --------------------------------------------------------------------

bim <- read.table(fileBim, header = FALSE, stringsAsFactors = FALSE)

# make unique names
bim$V2 <- 
  make.unique(
    ifelse(bim$V2 == ".",
           paste(bim$V1, bim$V4, bim$V5, bim$V6, sep = "_"),
           bim$V2),
    sep = "_")
            
# Output ------------------------------------------------------------------
write.table(bim, file = fileBim,
            col.names = FALSE, row.names = FALSE, quote = FALSE,
            sep = "\t")
