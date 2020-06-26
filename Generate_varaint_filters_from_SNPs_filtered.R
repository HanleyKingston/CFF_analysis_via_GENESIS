library(SeqArray)
library(SeqVarTools)

#snp_filtered_bi.tsv is only bi-allelic SNPs that pass GATK and VQSR filters
SNPs <- read.table("snp_filtered_bi.tsv", sep = "\t", header = TRUE)                 

nrow(SNPs)
#[1] 84360856

#spot_check that variants match GDS file by position and chromosome - must be TRUE!
gds <- seqOpen("CFF_sid_onlyGT.gds")

variants_gds <- data.frame(
    position = seqGetData(gds, "position"),
    chromosome = seqGetData(gds, "chromosome"),
    variant_id_gds = seqGetData(gds, "variant.id"),
    stringsAsFactors = FALSE
  )

#Run this multiple times
rand <- sample(1:nrow(SNPs),1)
SNPs[SNPs$variant.id == rand, "pos"] == variants_gds[variants_gds$variant_id_gds == rand, "position"] & SNPs[SNPs$variant.id == rand, "chrom"] == variants_gds[variants_gds$variant_id_gds == rand, "chromosome"] 


#Get a list of variants to use in LD-Pruning
saveRDS(SNPs$variant.id, "SNPS_bi_GATK_VQSR.rds")


#et a list of variants to use in association testing
#For association testing, also need to filter by missingness and MAF (this is given as an input to LD-pruning)
##Optional: plot minor allele frequency - note: PC-pruning also has an option to filter by MAF, so I will not include this in the filter
library(SeqVarTools)
gds <- seqOpen("CFF_sid_onlyGT.gds")


seqResetFilter(gds) #If there is a filter on the gds file, this will fail spectacularly

afreq <- alleleFrequency(gds) 
maf <- pmin(afreq, 1-afreq)
sum(maf <= 0.05)
#[1] 110327959
sum(maf > 0.05)
#[1] 9811885 #To keep in association testing


miss <- missingGenotypeRate(gds)
sum(miss < 0.05)
#To keep in association testing



maf_keep <- seqGetData(gds, "variant.id")[maf > 0.05]



#To test filter:
SNPS_bi_GATK_VQSR <- readRDS("SNPS_bi_GATK_VQSR.rds")
seqSetFilter(gds, variant.id = SNPS_bi_GATK_VQSR)

SNPS_bi_GATK_VQSR_MAF0.05_Miss0.05 <- readRDS("SNPS_bi_GATK_VQSR_MAF0.05_Miss0.05.rds")
seqSetFilter(gds, variant.id = SNPS_bi_GATK_VQSR_MAF0.05_Miss0.05)
