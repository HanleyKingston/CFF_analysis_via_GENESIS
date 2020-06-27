# LD-pruning

library(SeqArray)
library(SNPRelate)

# open GDS file
gds <- seqOpen("CFF_sid_onlyGT.gds") #This is the merged gds file ("CFF_5134_onlyGT.gds") from the seqArray_onlyGT folder, but I changed the sample.id's from the VCF_ids to the sids
variant_id <- readRDS("pre_LD_SNP_filter.rds")
#QC: https://github.com/HanleyKingston/CFF_analysis_via_GENESIS/blob/master/Generate_varaint_filters_from_SNPs_filtered.R
#Long range LD regions and chr7 excluded: https://github.com/HanleyKingston/CFF_analysis_via_GENESIS/blob/master/exclude_regions_beforeLD_prune.R
sample_id <- readRDS("keep_samples.rds")
#https://github.com/HanleyKingston/CFF_analysis_via_GENESIS/blob/master/generate_sample_filter_and_phenotype_df.R

# run LD pruning
snpset <- snpgdsLDpruning(gds,
                          sample.id = sample_id,
                          snp.id = variant_id,
                          maf = 0.05,
                          missing.rate = 0.05,
                          method = "corr",
                          slide.max.bp = 1 * 1e6, 
                          ld.threshold = sqrt(0.1)
                          autosome.only - TRUE
                          )


# convert list with one element per chrom to vector
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)
#[1] 171975

seqResetFilter(gds)
## of selected samples: 5,134
## of selected variants: 120,139,844

saveRDS(pruned, "6_26_prunedSNPs.rds")

# King (restart R)
library(SeqArray)
library(SNPRelate)

gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("6_26_prunedSNPs.rds")
sample_id <- readRDS("keep_samples.rds")

king <- snpgdsIBDKING(gds, snp.id = variant_id, sample.id = sample_id,
                      type = "KING-robust", autosome.only = TRUE)

kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id

saveRDS(king, paste0("6_26", "king_obj.rds"))
saveRDS(kingMat, paste0("6_26", "king_grm.rds")) #Note: should save a sparse matrix if using for association testing but not to give to PC-AiR



# Plot King (restart R)
library(ggplot2)
library(SNPRelate)

rel <- readRDS("6_26king_obj.rds")

kinship.df <- snpgdsIBDSelection(rel)

png(paste0("6_26_king", "_kinship.png"))
ggplot(kinship.df, aes_string("IBS0", "kinship")) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) + #Note: if you get a "partial transparancy is not supported..." error remove "alpha" argument
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()

100*sum(rel$kinship < 0)/length(rel$kinship)
#[1] 86.80962



# 1st itertion PC-Air (restart R)
library(SeqArray)
library(GENESIS)

gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("6_26_prunedSNPs.rds")
sample_id <- readRDS("keep_samples.rds")
kingMat <- readRDS("6_26king_grm.rds")

#NOTE: pcairpartition results are from an earlier run (actual varaints in LD_pruned will differ slightly due to randomness)
pc_part1 <- pcairPartition(gds, kinobj = kingMat, kin.thresh = 2^(-3), div.thresh = -2^(-4.5), divobj = kingMat)

str(pc_part1)
# $ rels  : chr [1:952] "S48216" "S71706" "S95760" "S11026" ...
# $ unrels: chr [1:4019] "S66366" "S48371" "S30713" "S43109" ...

mypcair <- pcair(gds, kinobj = kingMat, kin.thresh = 2^(-4.5),
                 divobj = kingMat, snp.include = variant_id,
                 sample.include = sample_id, div.thresh = -2^(-4.5))
#Working space: 3,979 samples, 171,975 SNVs #Unrelateds


saveRDS(mypcair, paste0("6_26_1it", "pcair.rds"))


# plot 1st iteration PCs with script
Rscript pca_plots.R 6_26_1itpcair.rds --out_prefix 6_26_1it --phenotype_file annot.rds --group race_or_ethnicity

#Plot correaltion by PCs
qsub -q new.q calculate_snp_pc_corr.sh
qsub -q new.q plot_snp_pc_corr.sh #Not working yet... need to install Topmed Pipeline or SSH into Adrienne's folder



# 1st ititeration PC-Relate (restart R)
library(SeqArray)
library(GENESIS)
library(SeqVarTools)

gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("6_26_prunedSNPs.rds")
sample_id <- readRDS("keep_samples.rds")
mypcair <- readRDS("6_26_1itpcair.rds")

seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id)

seqData <- SeqVarData(gds)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
str(iterator)

mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, seq(4)],
                    training.set = mypcair$unrels, sample.include = sample_id)



saveRDS(mypcrel, paste0("6_26_1it", "pcr_obj.rds"))
pcr_mat <- pcrelateToMatrix(mypcrel, scaleKin = 1)

saveRDS(pcr_mat, paste0("6,18_1it", "pcr_mat.rds"))



# Plot 1st ititeration PC-Relate (restart R)
library(ggplot2)
library(SNPRelate)

rel <- readRDS("6_26_1itpcr_obj.rds")

kinship.df <- rel$kinBtwn

png(paste0("6_26", "1it_kinship.png"))
ggplot(kinship.df, aes_string("k0", "kin")) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) + #Note: if you get a "partial transparancy is not supported..." error remove "alpha" argument
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()



# 2nd itertion PC-Air (restart R)
library(SeqArray)
library(GENESIS)

gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("pruned_excludedRegions_andChr7.rds")
sample_id <- readRDS("keep_samples.rds")
kingMat <- readRDS("6_26king_grm.rds")
pcrelate_matrix <- readRDS("6_26_1itpcr_mat.rds")

pc_part1 <- pcairPartition(gds, kinobj = pcrelate_matrix, kin.thresh = 2^(-4.5), div.thresh = -2^(-4.5), divobj = kingMat)
str(pc_part1)


mypcair <- pcair(gds, kinobj = pcrelate_matrix, kin.thresh = 2^(-4.5),
                 divobj = kingMat, snp.include = variant_id,
                 sample.include = sample_id, div.thresh = -2^(-4.5))

saveRDS(mypcair, paste0("6_26", "pcair.rds"))



# plot 1st iteration PCs with script
Rscript pca_plots.R 6_26pcair.rds --out_prefix 6_26 --phenotype_file annot.rds --group race_or_ethnicity

#Plot correaltion by PCs
qsub -q new.q calculate_snp_pc_corr.sh
qsub -q new.q plot_snp_pc_corr.sh #Not working yet... need to install Topmed Pipeline or SSH into Adrienne's folder




# 2nd ititeration PC-Relate
library(SeqArray)
library(GENESIS)
library(SeqVarTools)

gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("pruned_excludedRegions_andChr7.rds")
sample_id <- readRDS("keep_samples.rds")
mypcair <- readRDS("6_26_1itpcair.rds")

seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id)
## of selected samples: 4,971
## of selected variants: 171,679
seqData <- SeqVarData(gds)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)


mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, seq(4)],
                    training.set = mypcair$unrels, sample.include = sample_id)


saveRDS(mypcrel, paste0("6_26", "pcr_obj.rds"))
pcr_mat <- pcrelateToMatrix(mypcrel, scaleKin = 2)

saveRDS(pcr_mat, paste0("6_26", "pcr_mat.rds"))



# Plot 2nd ititeration PC-Relate
library(ggplot2)
library(SNPRelate)

rel <- readRDS("6,18pcr_obj.rds")

kinship.df <- rel$kinBtwn

png(paste0("6_26", "kinship.png"))
ggplot(kinship.df, aes_string("k0", "kin")) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) + #Note: if you get a "partial transparancy is not supported..." error remove "alpha" argument
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()

