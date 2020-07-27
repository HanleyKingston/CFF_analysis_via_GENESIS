library(SeqArray)
library(SeqVarTools)

#snp_filtered_bi.tsv is only bi-allelic SNPs that pass GATK and VQSR filters
SNPs <- read.table("/labdata12/CF_WGS2/shared/variants/snp_filtered_bi.tsv", sep = "\t", header = TRUE)                 

nrow(SNPs)
#[1] 87106997

#spot_check that variants match GDS file by position and chromosome - must be TRUE!
gds <- seqOpen("/labdata12/CF_WGS2/shared/variants/CFF_5134_onlyGT.gds")

variants_gds <- variantInfo(gds)



#Run this multiple times to verify gds variant IDs and snp_fitlered_bi variant IDs match by chromosome and position
rand <- sample(1:nrow(SNPs),1)
SNPs[SNPs$variant.id == rand, "pos"] == variants_gds[variants_gds$variant.id == rand, "pos"] & SNPs[SNPs$variant.id == rand, "chrom"] == variants_gds[variants_gds$variant.id == rand, "chr"] 

SNPs_for_LD_pruning <- SNPs$variant.id

#Get a list of variants to use in LD-Pruning
saveRDS(SNPs_for_LD_pruning, "SNPS_bi_GATK_VQSR.rds")


#Get a list of variants to use in association testing
#For association testing, also need to filter by missingness and MAF (this is given as an input to LD-pruning)
##Optional: plot minor allele frequency - note: PC-pruning also has an option to filter by MAF, so I will not include this in the filter

seqResetFilter(gds) #Note: if there is a filter on the gds file, this will fail spectacularly

afreq <- alleleFrequency(gds) 
maf <- pmin(afreq, 1-afreq)
sum(maf <= 0.05)
#[1] 110327959
sum(maf > 0.05)
#[1] 9811885 #To keep in association testing


miss <- missingGenotypeRate(gds)
sum(miss < 0.05)
#[1] 114432954 #To keep in association testing



maf_keep <- seqGetData(gds, "variant.id")[maf > 0.05]
miss_keep <- seqGetData(gds, "variant.id")[miss < 0.05]

SNPs_for_assoc_test <- intersect(SNPs_for_LD_pruning, intersect(maf_keep, miss_keep))


length(SNPs_for_assoc_test)
#[1] 5658280


saveRDS(SNPs_for_assoc_test, "SNPS_bi_GATK_VQSR_MAF0.05_miss0.05.rds")



#To test filter:
SNPS_bi_GATK_VQSR <- readRDS("SNPS_bi_GATK_VQSR.rds")
seqSetFilter(gds, variant.id = SNPS_bi_GATK_VQSR)
## of selected variants: 87,106,997

SNPS_bi_GATK_VQSR_MAF0.05_miss0.05 <- readRDS("SNPS_bi_GATK_VQSR_MAF0.05_miss0.05.rds")
seqSetFilter(gds, variant.id = SNPS_bi_GATK_VQSR_MAF0.05_miss0.05)
## of selected variants: 5,658,280


#get MAF ( >0.01) and MISS (<0.02) to match Paul's
seqResetFilter(gds) #Note: if there is a filter on the gds file, this will fail spectacularly


sum(maf > 0.01)
#[1] 14841927 #To keep in association testing


miss <- missingGenotypeRate(gds)
sum(miss < 0.02)
#[1] 113288995 #To keep in association testing



maf_keep <- seqGetData(gds, "variant.id")[maf > 0.01]
miss_keep <- seqGetData(gds, "variant.id")[miss < 0.02]

SNPs_for_assoc_test_match_Paul <- intersect(SNPs_for_LD_pruning, intersect(maf_keep, miss_keep))


length(SNPs_for_assoc_test_match_Paul)
#[1] 8186172


saveRDS(SNPs_for_assoc_test_match_Paul, "SNPS_bi_GATK_VQSR_MAF0.01_miss0.02.rds")

