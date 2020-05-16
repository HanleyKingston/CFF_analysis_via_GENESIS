#This generates a vector of sample IDs to be used as a filter in LD pruning, PC and GRM generation, and association testing 

library(SeqVarTools)

gds <- seqOpen("/home/hkings/DATA/CFF_5134_onlyGT.gds")
gds.id <- seqGetData(gds, "sample.id")

phen <- read.table("phenotype.txt", header = TRUE)

sum(phen$vcf_id %in% gds.id)
#[1] 5134
sum(gds.id %in% phen$vcf_id)
#[1] 5134
#Both columns must have all the same smaple IDs
#Check that they are also in the same order... THIS MUST BE TRUE
identical(as.character(phen$vcf_id), as.character(gds.id))

#Mine were not in the same order!, so must do this (can only do if the %in% checks are both equal):
phen <- phen[match(gds.id, phen$vcf_id),]
identical(as.character(phen$vcf_id), as.character(gds.id))


#Get a boolean vector of samples to exclude based on QC in phenotype data
#Must make sure to exclude all individuals who whould be excluded based on phenotype QC
sample_QC <- !is.na(phen$include_in_analysis)
table(sample_QC)
#FALSE  TRUE
#  120  5014


#To reset filter on gds, use
seqResetFilter(gds)

sample_filter <- seqGetData(gds, "sample.id")[sample_QC]

saveRDS(sample_filter, file = "keep_samples.rds")

