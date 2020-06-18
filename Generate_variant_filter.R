#This generates a list of variant IDs based on filter criteria... I generated a more and less stringent filter
#both filters follows suggested filters from: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112 (see:[C] Hard-filter SNPs on multiple expressions using VariantFiltration)
#My gds file already excldues the X chromosome
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


flag.metric.df <- readRDS("flag.metric_with_gds_ids.rds")

#alt:
#load("flag.metric.RData")

str(flag.metric.df)


#filter criteria:

table(flag.metric.df$flag.info_QD.2, useNA = "ifany")
#    FALSE      TRUE
#114175552   5964292
QD <- flag.metric.df$flag.info_QD.2 == FALSE

table(flag.metric.df$flag.qual.30)
#    FALSE
#120139844
qual <- flag.metric.df$flag.qual.30 == FALSE

table(flag.metric.df$flag.info_SOR.3)
#    FALSE      TRUE
#116592646   3547198
SOR <- flag.metric.df$flag.info_SOR.3 == FALSE

table(flag.metric.df$flag.info_FS.60)
#    FALSE      TRUE
#118469110   1670734
FS <- flag.metric.df$flag.info_FS.60 == FALSE

table(flag.metric.df$flag.info_MQ.40)
#    FALSE      TRUE
#110211257   9928587
MQ <- flag.metric.df$flag.info_MQ.40 == FALSE

table(flag.metric.df$flag.info_ReadPosRankSum.8)
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

table(flag.metric.df$flag.missByVar.0.0, useNA = "ifany"5)
#    FALSE      TRUE
#114437649   5702195
missByVar <- flag.metric.df$flag.missByVar.0.05 == FALSE

#Exclude X-chromosome (this is probably happens based on other filters, but just to be safe)
table(flag.metric.df$chr == "X")
#    FALSE      TRUE
#115537943   4601901
noX <- flag.metric.df$chr != "X"


##Get minor allele frequency - note: PC-pruning also has an option to filter by MAF, but I am doing it here for consistency with variants I will use in association testing)
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
sum(flag.metric.df$variant_id_gds) == nrow(flag.metric.df$variant_id_gds) #Check that all varaint IDs are unique
var_filter_moderate_gds.temp <- flag.metric.df[flag.metric.df$variant.id, "variant_id_gds"]


## Extract the correct varaint ID... IMPORTANT: must make sure gds varaint IDs and flag.metric IDs match by position and chromosome (see "Check_or_match_gdsIDs_to_Filters.R")

### Moderate
var_filter_moderate_gds.temp <- flag.metric.df$variant_id_gds[RankSum & MQ & FS & SOR & qual & QD & (snvPASS | indelPASS) & noX] 
var_filter_SNPs_and_indels_MAF0.01 <- intersect(seqGetData(gds, "variant.id")[maf > 0.01], var_filter_moderate_gds.temp)
length(var_filter_SNPs_and_indels_MAF0.01)
#[1] 10418664

### Stringent
var_filter_stringent_gds.temp <- flag.metric.df$variant_id_gds[RankSum & MQ & FS & SOR & qual & QD & snvPASS & missByVar & biAllelic & noX]
var_filter_SNVs_MAF0.05 <- intersect(seqGetData(gds, "variant.id")[maf > 0.05], var_filter_stringent_gds.temp)
length(var_filter_SNVs_MAF0.05)
#[1] 5490945


### Just VQSR & GATK hard filters:
VQSR_GATK_only <- flag.metric.df$variant_id_gds[RankSum & MQ & FS & SOR & qual & QD & (snvPASS | indelPASS)] 


#Note: not filtering by HW (should have already been done and I don't have the family info to control for relatedness yet)


saveRDS(var_filter_SNPs_and_indels_MAF0.01, file = "var_filter__SNVs_and_indels_lowMAF.rds")
saveRDS(var_filter_SNVs_MAF0.05, file = "var_filter_SNVs_MAF0.05.rds")
saveRDS(VQSR_GATK_only, file = "var_filter_VQSR_GATK_only.rds")


#To test filter:
seqSetFilter(gds, variant.id = var_filter_SNVs_MAF0.05)
## of selected variants: 5,490,945






