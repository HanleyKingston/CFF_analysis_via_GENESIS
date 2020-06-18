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


##Get minor allele frequency
library(SeqVarTools)

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


flag.metric.df <- readRDS("flag.metric_with_gds_ids.rds")

#alt:
#load("flag.metric.RData")

str(flag.metric.df)




#filter criteria:

table(flag.metric.df$flag.info_QD.2, useNA = TRUE)
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


gds <- seqOpen("/home/hkings/DATA/CFF_sid_onlyGT.gds")

length(seqGetData(gds, "variant.id"))
#[1] 115537943
nrow(flag.metric.df)
#[1] 120139844
sum(seqGetData(gds, "variant.id") %in% flag.metric.df$variant.id)
#[1] 115537943
sum(seqGetData(gds, "variant.id") %in% gds.var)
#[1] 115537943
#The difference is probably because I excluded the X chromosome from my merged gds file


#Make a datframe of gds variant info
position <- seqGetData(gds, "position")
chromosome <- seqGetData(gds, "chromosome")
variant <- seqGetData(gds, "variant.id")
gds.ids <- data.frame(position, chromosome, variant, stringsAsFactors=FALSE)

head(gds.ids)

#Exclude X chromosome from flag.metric.df because not in my gds file
flag.metrix.df2 <- flag.metric.df[flag.metric.df$chr != "X",]

#dataframes are now the same size
nrow(gds.ids)
#[1] 115537943
nrow(flag.metric.df2)
#[1] 115537943

#Add a column for gds variant ID to flag.metric dataframe
flag.metrix.df2$gds.var.id <- NULL 

for(i in 1:nrow(flag.metric.df2)){
    flag.metric.df2$gds.var.id[i] <- gds.ids[gds.ids$chromosome == flag.metric.df2$chr[i] & gds.ids$position == flag.metric.df2$pos[i], "variant"]
    if(i == 1000000){
        print("1 mil run")
      }
}

gds.ids[

#This didn't work
head(flag.metric.df2)
#...
#  flag.info_RAWMQ.50000 flag.info_RAWMQ.1e+05 gds.var.id
#1                 FALSE                 FALSE          1
#2                 FALSE                 FALSE          2
#3                 FALSE                 FALSE          3

#The gds.var.id is just 1 after the 1st chromosome:
gdsflag.metric.df2[match("2",flag.metric.df2$chr),]
#...
#        flag.info_RAWMQ.50000 flag.info_RAWMQ.1e+05 gds.var.id
#39936670                 FALSE                 FALSE          1
flag.metric.df2[39936670:39936675,]
#         gds.var.id
#39936670          1
#39936671          1
#39936672          1

gds.ids[match("2", gds.ids$chromosome),]
#...
#         position chromosome  variant
#51629765    10188          2 51629765

flag.metric.df2[match("3",flag.metric.df2$chr),]
#...
#         gds.var.id
#19814256          1


for(i in unique(gds.ids$chromosome)){
  gds.ids.chr <- gds.ids[gds.ids$chromosome == i,]
  print(gds.ids[duplicated(gds.ids.chr$position),])
  }

gds.ids[gds.ids$chromosome == 
        

gds.ids <- gds.ids[order(as.numeric(chromosome)),]
gds.ids[match("2", gds.ids$chromosome),]
#         position chromosome  variant
#51629765    10188          2 51629765
gds.ids[51629764,]

        


library(stringr)
gds.ids <- gds.ids[str_sort(gds.ids$chromosome, numeric = TRUE),]
gds.ids[match("2", gds.ids$chromosome),]


for(i in 1:nrow(gds.ids)){
    flag.metric.df2$gds.var.id[i] <- gds.ids[gds.ids$chromosome == flag.metric.df2$chr[i] & gds.ids$position == flag.metric.df2$pos[i], "variant"]
}

var_filter_moderate.temp <- flag.metric.df$gds.var.id[RankSum & MQ & FS & SOR & qual & QD & (snvPASS | indelPASS)] 
#Convert to gds varaints if necissary
sum(flag.metric.df$variant.id) == nrow(flag.metric.df$variant.id) #Check that all varaint IDs are unique
sum(flag.metric.df$variant_id_gds) == nrow(flag.metric.df$variant_id_gds) #Check that all varaint IDs are unique
var_filter_moderate_gds.temp <- flag.metric.df[flag.metric.df$variant.id, "variant_id_gds"]
#Include MAF for final filter
var_filter_SNPs_and_indels_MAF0.01 <- intersect(seqGetData(gds, "variant.id")[maf > 0.01], var_filter_moderate_gds.temp)
length(var_filter__SNVs_and_indels_lowMAF)
#[1] 11864290

var_filter_stringent.temp <- flag.metric.df$gds.var.id[RankSum & MQ & FS & SOR & qual & QD & snvPASS & missByVar & biAllelic]
#Convert to gds varaints if necissary
var_filter_moderate_gds.temp <- flag.metric.df[flag.metric.df$variant.id, "variant_id_gds"]
var_filter_SNVs_MAF0.5 <- intersect(seqGetData(gds, "variant.id")[maf > 0.05], var_filter_stringent_gds.temp)
length(var_filter_SNVs_MAF0.05)
#[1] 6829196

#Note: not filtering by HW (should have already been done and I don't have the family info to control for relatedness yet)


saveRDS(var_filter_moderate, file = "var_filter__SNVs_and_indels_lowMAF.rds")
saveRDS(var_filter_stringent, file = "var_filter_SNVs_MAF0.05")

#To test filter:
seqSetFilter(gds, variant.id = var_filter_stringent)
## of selected variants: 6,829,196





