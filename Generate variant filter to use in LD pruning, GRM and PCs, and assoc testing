#The genrates a list of variant IDs based on filter criteria... I generated a more and less stringent filter
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


library(SeqVarTools)

gds <- seqOpen("/home/hkings/DATA/CFF_5134_onlyGT.gds")

#To reset filter on gds, use
seqResetFilter(gds)

##Get minor allele frequency
afreq <- alleleFrequency(gds) #Can also get this from hw$afreq
maf <- pmin(afreq, 1-afreq)
sum(maf <= 0.05)
#[1] 106032277 
sum(maf > 0.05)
#[1] 9505666 #These I will keep

#Plot MAF
pdf("minor_allele_freq_hist.pdf")
hist(maf, breaks=50)
dev.off()

load("flag.metric.RData")
str(flag.metric.df)



#filter criteria:

table(flag.metric.df$flag.info_QD.2)
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

table(flag.metric.df$index.snvPASS)
#   FALSE     TRUE
#27543319 92
6525
snvPASS <- flag.metric.df$index.snvPASS == TRUE

table(flag.metric.df$index.indelPASS)
#    FALSE      TRUE
#107380642  12759202
indelPASS <- flag.metric.df$index.indelPASS == TRUE

table(flag.metric.df$index.snv_biAllelic)
#    FALSE      TRUE
# 19202784 100937060
biAllelic <- flag.metric.df$index.snv_biAllelic == TRUE

table(flag.metric.df$flag.missByVar.0.05)
#    FALSE      TRUE
#114437649   5702195
missByVar <- flag.metric.df$flag.missByVar.0.05 == FALSE


gds.var <- seqGetData(gds, "variant.id")
length(gds.var)
#[1] 115537943
nrow(flag.metric.df)
[1] 120139844
sum(gds.var %in% flag.metric.df$variant.id)
#[1] 115537943
sum(flag.metric.df$variant.id %in% gds.var)
#[1] 115537943
#The difference is probably because I excluded the X chromosome from my merged gds file


var_filter_moderate.temp <- flag.metric.df$variant.id[RankSum & MQ & FS & SOR & qual & QD & (snvPASS | indelPASS)]
var_filter_moderate <- intersect(seqGetData(gds, "variant.id")[maf > 0.01], var_filter_moderate.temp)
length(var_filter_moderate)
#[1] 11864290

var_filter_stringent.temp <- flag.metric.df$variant.id[RankSum & MQ & FS & SOR & qual & QD & snvPASS & missByVar & biAllelic]
var_filter_stringent <- intersect(seqGetData(gds, "variant.id")[maf > 0.05], var_filter_stringent.temp)
length(var_filter_stringent)
#[1] 6829196

#Note: not filtering by HW (should have already been done and I don't have the family info to control for relatedness yet)


saveRDS(var_filter_moderate, file = "keep_var_moderate.rds")
saveRDS(var_filter_stringent, file = "keep_var_stringent.rds")

#To test filter:
seqSetFilter(gds, variant.id = var_filter_stringent)





