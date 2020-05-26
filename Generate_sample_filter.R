#This generates a vector of sample IDs to be used as a filter in LD pruning, PC and GRM generation, and association testing 

library(SeqVarTools)

gds <- seqOpen("/home/hkings/DATA/CFF_sid_onlyGT.gds")
gds.id <- seqGetData(gds, "sample.id")

phenotype <- read.table("phenotype.txt", header = TRUE)

##Both data frames must have all the same smaple IDs
sum(phen$sid %in% gds.id)
#[1] 5134
sum(gds.id %in% phen$sid)
#[1] 5134


#THIS MUST BE TRUE
identical(as.character(phenotype$sid), as.character(gds.id))

#These should be in the same order from the last script, but if not, must do this (can only do if the %in% checks are both equal):
#phenotype <- phenotype[match(gds.id, phenotype$sid),]

#Get a boolean vector of samples to exclude based on QC in phenotype data
#Must make sure to exclude all individuals who whould be excluded based on phenotype QC
sample_QC <- !is.na(phenotype$include_in_analysis)
table(sample_QC)
#FALSE  TRUE
#  165  4971


#To reset filter on gds, use
seqResetFilter(gds)

sample_filter <- seqGetData(gds, "sample.id")[sample_QC]

saveRDS(sample_filter, file = "keep_samples.rds")

