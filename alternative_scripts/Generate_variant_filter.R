#This generates a list of variant IDs based on filter criteria... I generated a more and less stringent filter
#both filters follows VQSR & GATK hard filters:suggested filters from: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112 (see:[C] Hard-filter SNPs on multiple expressions using VariantFiltration)
#
#moderate filter also excludes:
# MAF < 0.01
# SNPs and indels that don't pass VQSR filters (index.snvPASS = FALSE | index.indelPASS = FALSE)
#
#
#stringent filter also excludes:
# MAF < 0.05
# missingness by variant > 0.05
# indels
# anything not biallelic
# SNPs that don't pass VQSR filters (index.snvPASS = FALSE)
#
#Note: can also filter by MAF and missingness in GENESIS's LD-pruning, but I chose to do it here so I can use the same filter for the association testing


flag.metric.df <- readRDS("/home/hkings/DATA/old_files/flag.metric_with_gds_ids.rds")

#alt:
#load("flag.metric.RData")

str(flag.metric.df)


#filter criteria:
#Note: Suggested filter, but not available: MQRankSum > 12.5

table(flag.metric.df$flag.info_QD.2, useNA = "ifany")
#    FALSE      TRUE
#114175552   5964292
QD <- flag.metric.df$flag.info_QD.2 == FALSE

table(flag.metric.df$flag.qual.30, useNA = "ifany")
#    FALSE
#120139844
qual <- flag.metric.df$flag.qual.30 == FALSE

table(flag.metric.df$flag.info_SOR.3, useNA = "ifany")
#    FALSE      TRUE
#116592646   3547198
SOR <- flag.metric.df$flag.info_SOR.3 == FALSE

table(flag.metric.df$flag.info_FS.60, useNA = "ifany")
#    FALSE      TRUE
#118469110   1670734
FS <- flag.metric.df$flag.info_FS.60 == FALSE

table(flag.metric.df$flag.info_MQ.40, useNA = "ifany")
#    FALSE      TRUE
#110211257   9928587
MQ <- flag.metric.df$flag.info_MQ.40 == FALSE

table(flag.metric.df$flag.info_ReadPosRankSum.8, useNA = "ifany")
#    FALSE      TRUE
#120137053      2791
RankSum <- flag.metric.df$flag.info_ReadPosRankSum.8 == FALSE

table(flag.metric.df$index.snvPASS, useNA = "ifany")
#   FALSE     TRUE
#27543319 92659625
snvPASS <- flag.metric.df$index.snvPASS == TRUE

table(flag.metric.df$index.indelPASS, useNA = "ifany")
#    FALSE      TRUE
#107380642  12759202
indelPASS <- flag.metric.df$index.indelPASS == TRUE

table(flag.metric.df$index.snv_biAllelic, useNA = "ifany")
#    FALSE      TRUE
# 19202784 100937060
biAllelic <- flag.metric.df$index.snv_biAllelic == TRUE

table(flag.metric.df$flag.missByVar.0.5, useNA = "ifany")
#    FALSE      TRUE
#114437649   5702195
missByVar <- flag.metric.df$flag.missByVar.0.05 == FALSE

#Exclude X-chromosome
table(flag.metric.df$chr == "X")
#    FALSE      TRUE
#115537943   4601901
noX <- flag.metric.df$chr != "X"


##Optional: plot minor allele frequency - note: PC-pruning also has an option to filter by MAF, so I will not include this in the filter
library(SeqVarTools)
gds <- seqOpen("CFF_sid_onlyGT.gds")

afreq <- alleleFrequency(gds) #Can also get this from hw$afreq
maf <- pmin(afreq, 1-afreq)
sum(maf <= 0.05)
#[1] 110327959
sum(maf > 0.05)
#[1] 9811885 #These I will keep

#Plot MAF
png("minor_allele_freq_hist.png")
hist(maf, breaks=50)
dev.off()


length(seqGetData(gds, "variant.id"))
#[1] 120139844
nrow(flag.metric.df)
#[1] 120139844

#If file has multiple varaint ID's, make sure there are no repeats:
length(unique(flag.metric.df$variant_id)) == length(flag.metric.df$variant_id) #Check that all varaint IDs are unique


## Extract the correct varaint ID... IMPORTANT: must make sure gds varaint IDs and flag.metric IDs match by position and chromosome (see "Check_or_match_gdsIDs_to_Filters.R")

### Moderate
var_filter_SNVs_and_indels <- flag.metric.df$variant_id[RankSum & MQ & FS & SOR & qual & QD & (snvPASS | indelPASS)] #Add in missingness
length(var_filter_SNVs_and_indels)
#[1] 100635332


### Stringent
var_filter_SNVs <- flag.metric.df$variant_id[RankSum & MQ & FS & SOR & qual & QD & snvPASS & missByVar & biAllelic]
length(var_filter_SNVs)
#[1] 87672807


#Note: not filtering by HW (should have already been done and I don't have the family info to control for relatedness yet)


saveRDS(var_filter_SNVs_and_indels , file = "var_filter_SNVs_and_indels.rds")
saveRDS(var_filter_SNVs, file = "var_filter_SNVs.rds")


#To test filter:
seqSetFilter(gds, variant.id = var_filter_SNVs)
## of selected variants: 5,490,945
seqResetFilter(gds)
## of selected variants: 87,672,807






