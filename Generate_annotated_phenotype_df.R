
#Load phenotype data:
phenotype <- read.table("phenotype.txt", sep = "\t", header = TRUE)

##The null model requires a sample ID column:
phenotype$sample.id <- as.character(phenotype$sid)

##Check that columns match gds file
library(SeqVarTools)
gdsfmt::showfile.gds(closeall=TRUE)
gds <- seqOpen("CFF_sid_onlyGT.gds")
gds.id <- seqGetData(gds, "sample.id")
sum(gds.id %in% phenotype$sample.id) == sum(phenotype$sample.id %in% gds.id) #phenotype file isn't fitlered so will have more total samples
keep_samples <- readRDS("keep_samples.rds")


identical(as.character(phenotype[phenotype$sample.id %in% keep_samples, "sample.id"]), gds.id)
