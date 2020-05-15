# CFF_analysis_via_GENESIS

#To access R:
module load conda
R --version
R

#To access python:
module load conda/5.3.1
source activate py3
python --version
python


###R###
#Save gds sample.id as a file
library(SeqArray)

gdsfile <- "/labdata12/CF_WGS2/Variants/CFF_5134_GDSs/seqArray_onlyGT/CFF_5134_chr21_onlyGT.gds"
gds <- seqOpen(gdsfile)
gds.id <- seqGetData(gds, "sample.id")
length(gds.id)
#[1] 5134

head(gds.id)
[1] "JHU_BL_01_1" "JHU_BL_02_1" "JHU_BL_03_1" "JHU_BL_04_1" "JHU_BL_05_1" "JHU_BL_06_1"

write(gds.id, "/home/hkings/DATA/sample_id_gds.txt")

q()
### ###


###PYTHON###
import re

flags = open("samples_flag_or_exclude.txt", "r")

toss_flag = open("toss_flag_list.txt", "w")

for line in flags:
    match = re.search("(\| |\|)+(UNC|UW|JHU)+(\S*)", line)
    if match:
        if "Drop" in line or "?" in line:
            code = "Exclude"
            toss_flag.write(match.group(2)+match.group(3) + "\t" + code + "\n")
        elif "Keep" in line or "flagged" in line:
            code = "Flag"
            toss_flag.write(match.group(2)+match.group(3) + "\t" + code + "\n")
        else:
            toss_flag.write(match.group(2)+match.group(3) + "\n")
    else:
        toss_flag.write(line + "\n")

toss_flag.close()
### ###

#From key, I found one ind. that says hold_for_UW_ID (JHU_TSS_933_2641_1) I flagged this person by their sample #

#List of duplicates (by Participant_ID_corrected) to exclude
#cut -f4 key_cffwgs.tsv | sort | uniq -cd > non_unique_ID_corrected.txt


#Add Subject_ID_Corrected column to participants_cffwgs.tsv
#join -1 1 -2 2 -t $'\t' <(sort -k1,1 participants_cffwgs.tsv) <(sort -k2,2 key_cffwgs.tsv) > sorted.txt
#Edit to print out anyone who doesn't match between files (I already checked in R though that all pids are in both file... I think!

###R###
#subjects <- read.delim("sorted.txt", sep = "\t", header = TRUE)
participants <- read.delim("participants_cffwgs.tsv", sep = "\t", header = TRUE)
nrow(participants)
#[1] 5161
key <- read.delim("key_cffwgs.tsv", sep = "\t", header = TRUE)
nrow(key)
#[1] 5199
flag <- scan("samples_flag.txt", "character", sep = "\n")
#Read 79 items
exclude <- scan("samples_exclude.txt", "character", sep = "\n")
#Read 51 items
gds.id <- scan("sample_id_gds.txt", "character", sep = "\n")
#Read 5134 items

"%notin%" <- Negate("%in%")

sum(key$vcf_id %notin% gds.id)
#[1] 65
in_key_not_gds <- key[key$vcf_id %notin% gds.id,"vcf_id"]
#These should be removed!

sum(participants$pid %in% key$pid)
#[1] 5161 #All participants are in the key
sum(key$pid %in% participants$pid)
#[1] 5199 #All indiviudals in key are in participants, so 38 people must be repeated in the key
sum(duplicated(key$pid))
#[1] 38 #38 pairs of duplicates, so difference between key and participants file is just that duplicates were removed
#I only want to keep one of each duplicate... keep the person who is in the gds file
#Get all duplicates:
key[duplicated(key$pid) | duplicated(key$pid, fromLast = TRUE),]

#Because each duplicate has a different vcf_id and sid, I want to keep the one in the gds.id file
duplicates <- key[duplicated(key$pid) | duplicated(key$pid, fromLast = TRUE), "vcf_id"]
length(duplicates)
#[1] 72
sum(duplicates %in% gds.id)
#[1] 71
#all but one of the duplicates are in the gds file... since I can't determine which is the correct sample, I will remove them from the analysis entirely?
sum(in_key_not_gds %in% duplicates)
#[1] 1 #One of the IDs that are in the key but not the gds are in duplicates
sum(exclude %in% duplicates)
#[1] 0


#key[duplicated(key$pid) | duplicated(key$pid, fromLast = TRUE), c("sid", "vcf_id")] <- "NA"
nrow(key)
#[1] 5199
sum(key$vcf_id %notin% gds.id)
#[1] 65

#Add key and phenotype together so phenotype data contains the vcf_id
phenotype <- merge(participants, key, by = "pid", all = TRUE)

phenotype <- phenotype[order(phenotype$vcf_id),]
gds.id <- sort(gds.id)

sum(phenotype$vcf_id %notin% gds.id)
#[1] 65
not_shared_in_phenotype <- phenotype[phenotype$vcf_id %notin% gds.id,"vcf_id"]
sum(exclude %in% not_shared_in_phenotype)
#[1] 0
sum(flag %in% not_shared_in_phenotype)
#[1] 18


sum(gds.id %notin% phenotype$vcf_id)
#[1] 0 #There are no samples in the gds file without corresponding phenotype info

#anything in phenotype file without a corresponding entry in gds file should be removed, otherwise analysis won't run
phenotype_pruned <- phenotype[phenotype$vcf_id %in% gds.id,]
nrow(phenotype_pruned)
#[1] 5134
sum(phenotype_pruned$vcf_id %notin% gds.id)
#[1] 0


phenotype_pruned$include_in_analysis <- ifelse(phenotype_pruned$vcf_id %in% duplicates | phenotype_pruned$vcf_id %in% exclude, "NA", ifelse(phenotype_pruned$vcf_id %in% flag, "flag", "keep"))
table(phenotype_pruned$include_in_analysis)
#flag keep   NA
#  39 4975  120

colnames(phenotype_pruned)

#To make more managable, I'm just selecting phenotypes I'm interested in
phenotype_pruned_temp  <- phenotype_pruned[,c("pid", "sex_wgs", "birthdate_year", "cftr_var_1_wgs", "cftr_var_2_wgs", "cftr_addl_vars_wgs", "cftr_gt_category_wgs", "age_dx", "year_dx", "age_death", "knorma", "vcf_id", "include_in_analysis")]

write.table(phenotype_pruned, "phenotype.txt", sep = "\t")

### ###



Rscript mergeGDS_onlyGT.R #The way the script is set up, need to run this in R if the ind. chromosome files aren't in the same folder as where you want the merged file
#9,081,025,518 bytes in total
#5134 samples in total (5134 samples in common)
#115, 537,943 variants


###R### QC
library(SeqVarTools)

gds <- seqOpen("/home/hkings/DATA/CFF_5134_onlyGT.gds")
gds.id <- seqGetData(gds, "sample.id")

phen <- read.table("phenotype.txt", header = TRUE)

sum(phen$vcf_id %in% gds.id)
#[1] 5134
sum(gds.id %in% phen$vcf_id)
#[1] 5134
#Both columns have all the same smaple IDs!
identical(as.character(phen$vcf_id), as.character(gds.id))
#[1] FALSE

#Mine were not in the same order!, so must do this (can only do if the %in% checks are both equal):
phen <- phen[match(gds.id, phen$vcf_id),]
identical(as.character(phen$vcf_id), as.character(gds.id))
#[1] TRUE


#Get a boolean vector of samples to exclude based on QC in phenotype data
#Must make sure to exclude all individuals who whould be excluded based on phenotype QC
sample_QC <- !is.na(phen$include_in_analysis)
table(sample_QC)
#FALSE  TRUE
#  120  5014



#To reset filter on gds, use
seqResetFilter(gds)
# of selected samples: 5,134
# of selected variants: 115,537,943

sample_filter <- seqGetData(gds, "sample.id")[sample_QC]

saveRDS(sample_filter, file = "keep_samples.rds")




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
biAllelic <- flag.metric.df$index.indelPASS == FALSE

table(flag.metric.df$flag.missByVar.0.05)
#    FALSE      TRUE
#114437649   5702195
missByVar <- flag.metric.df$flag.missByVar.0.05 == FALSE




gds.var <- seqGetData(gds, "variant.id")]
length(gds.var)
#[1] 115537943

nrow(flag.metric.df)
[1] 120139844





sum(gds.var %in% flag.metric.df$variant.id)
#[1] 115537943
sum(flag.metric.df$variant.id %in% gds.var)
#[1] 115537943
#The difference is probably because I excluded the X chromosome from my merged gds file



var_filter_moderate.temp <- flag.metric.df$variant.id[RankSum & MQ & FS & SOR & qual & QD & missByVar & (snvPASS | indelPASS)]
var_filter_moderate <- intersect(seqGetData(gds, "variant.id")[maf > 0.01], var_filter_moderate.temp)
length(var_filter_moderate)
#[1] 11810424

var_filter_stringent.temp <- flag.metric.df$variant.id[RankSum & MQ & FS & SOR & qual & QD & snvPASS & missByVar & biAllelic]
var_filter_stringent <- intersect(seqGetData(gds, "variant.id")[maf > 0.05], var_filter_stringent.temp)
length(var_filter_stringent)
#[1] 7005599

#Note: not filtering by HW (should have already been done and I don't have the family info to control for relatedness yet)


saveRDS(var_filter_moderate, file = "keep_var_moderate.rds")
saveRDS(var_filter_stringent, file = "keep_var_stringent.rds")



seqSetFilter(gds, sample.id=sample_filter, variant.id = var_filter_stringent)
# of selected samples: 5,014
# of selected variants: 7,005,599





#Generate PCs and GRM using LD pruning threshold of sqrt(0.1)
Rscript PC_and_grm_script.R CFF_5134_onlyGT.gds 0.316227766 "keep_variants.rds" "keep_samples.rds" "LDsqrt0.1" & > LDsqrt0.1_PCs_grm_script.out




