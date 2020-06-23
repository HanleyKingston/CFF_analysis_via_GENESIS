# LD-pruning

library(SeqArray)
library(SNPRelate)

# open GDS file
gds <- seqOpen("CFF_sid_onlyGT.gds") #This is the merged gds file ("CFF_5134_onlyGT.gds") from the seqArray_onlyGT folder, but I changed the sample.id's from the VCF_ids to the sids
variant_id <- readRDS("var_filter_SNVs_MAF0.05.rds")
sample_id <- readRDS("keep_samples.rds")

# run LD pruning
snpset <- snpgdsLDpruning(gds,
                          sample.id = sample_id,
                          snp.id = variant_id,
                          maf = 0.05,
                          missing.rate = 0.05,
                          method = "corr",
                          slide.max.bp = 1 * 1e6, 
                          ld.threshold = sqrt(0.1)
                          )

#SNV pruning based on LD:
#Calculating allele counts/frequencies ...
#[==================================================] 100%, completed, 9.0m
#Excluding 6,828 SNVs (monomorphic: TRUE, MAF: 0.05, missing rate: 0.05)
#Working space: 4,971 samples, 5,484,117 SNVs
#    using 1 (CPU) core
#    sliding window: 1,000,000 basepairs, Inf SNPs
#    |LD| threshold: 0.316228
#    method: correlation
#Chromosome 1: 0.14%, 13,864/9,936,669
#Chromosome 10: 0.16%, 8,934/5,678,525
#Chromosome 11: 0.14%, 8,298/5,822,215
#Chromosome 12: 0.15%, 8,596/5,561,704
#Chromosome 13: 0.15%, 6,466/4,404,461
#Chromosome 14: 0.17%, 5,996/3,567,216
#Chromosome 15: 0.17%, 5,887/3,417,178
#Chromosome 16: 0.17%, 6,568/3,826,829
#Chromosome 17: 0.18%, 6,225/3,380,458
#Chromosome 18: 0.17%, 5,930/3,417,044
#Chromosome 19: 0.21%, 5,456/2,617,465
#Chromosome 2: 0.14%, 13,464/9,877,586
#Chromosome 20: 0.18%, 5,103/2,842,811
#Chromosome 21: 0.17%, 2,967/1,722,032
#Chromosome 22: 0.19%, 3,413/1,809,412
#Chromosome 3: 0.15%, 11,614/8,005,377
#Chromosome 4: 0.14%, 11,175/7,759,638
#Chromosome 5: 0.15%, 10,548/7,235,943
#Chromosome 6: 0.15%, 10,105/6,751,274
#Chromosome 7: 0.15%, 9,708/6,373,532
#Chromosome 8: 0.14%, 8,827/6,159,400
#Chromosome 9: 0.16%, 8,395/5,371,174
#177,539 markers are selected in total.


# convert list with one element per chrom to vector
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)
#[1] 177539

seqResetFilter(gds)
## of selected samples: 5,134
## of selected variants: 120,139,844

# Also filter out 4 problematic LD regions and excldue chromosome 7
filt <- get(data(list=paste("pcaSnpFilters", "hg38", sep="."), package="GWASTools"))

add_regions <- data.frame(t(c(7, 1, 159345973, "to exclude deltaF508")))
colnames(add_regions) <- colnames(filt)
rownames(add_regions) <- "chr7"
filt <- rbind(filt, add_regions)

chrom <- seqGetData(gds, "chromosome")
pos <- seqGetData(gds, "position")
pca.filt <- rep(TRUE, length(chrom))
for (f in 1:nrow(filt)) {
    pca.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos & pos < filt$end.base[f]] <- FALSE
}
seqSetFilter(gds, variant.sel=pca.filt, action="intersect", verbose=TRUE)
## of selected variants: 116,399,609
}
  
PCAcorr_snps <- seqGetData(gds, "variant.id")
saveRDS(PCAcorr_snps, "PCAcorr_snps.rds")

pruned_excludedRegions_andChr7 <- intersect(pruned, PCAcorr_snps)
length(pruned_excludedRegions_andChr7)
#[1] 171679

saveRDS(pruned_excludedRegions_andChr7, "pruned_excludedRegions_andChr7.rds")
write(pruned_excludedRegions_andChr7, "pruned_excludedRegions_andChr7.txt", ncol = 1) 



# King (restart R)
library(SeqArray)
library(SNPRelate)

gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("pruned_excludedRegions_andChr7.rds")
sample_id <- readRDS("keep_samples.rds")

king <- snpgdsIBDKING(gds, snp.id = variant_id, sample.id = sample_id,
                      type = "KING-robust")
#IBD analysis (KING method of moment) on genotypes:
#Calculating allele counts/frequencies ...
#[==================================================] 100%, completed, 6.8m
#Working space: 4,971 samples, 171,679 SNVs
#    using 1 (CPU) core
#No family is specified, and all individuals are treated as singletons.
#Relationship inference in the presence of population stratification.
#CPU capabilities: Double-Precision SSE2
#Thu Jun 18 15:58:59 2020    (internal increment: 19968)
#[==================================================] 100%, completed, 10.4m
#Thu Jun 18 16:09:25 2020    Done.


kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id

saveRDS(king, paste0("6,18", "king_obj.rds"))
saveRDS(kingMat, paste0("6,18", "king_grm.rds")) #Note: should save a sparse matrix if using for association testing but not to give to PC-AiR



# Plot King (restart R)
library(ggplot2)
library(SNPRelate)

rel <- readRDS("6,18king_obj.rds")

kinship.df <- snpgdsIBDSelection(rel)

png(paste0("6,18_king", "_kinship.png"))
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
variant_id <- readRDS("pruned_excludedRegions_andChr7.rds")
sample_id <- readRDS("keep_samples.rds")
kingMat <- readRDS("6,18king_grm.rds")

pc_part1 <- pcairPartition(gds, kinobj = kingMat, kin.thresh = 2^(-3), div.thresh = -2^(-4.5), divobj = kingMat)
#Using kinobj and divobj to partition samples into unrelated and related sets
#Identifying relatives for each sample using kinship threshold 0.125
#Identifying pairs of divergent samples using divergence threshold -0.0441941738241592
#Partitioning samples into unrelated and related sets...
#Warning message:
#In pcairPartition(gds, kinobj = kingMat, kin.thresh = 2^(-3), div.thresh = -2^(-4.5),  :
#  some samples in unrel.set are not in kinobj or divobj; they will not be included
str(pc_part1)
# $ rels  : chr [1:952] "S48216" "S71706" "S95760" "S11026" ...
# $ unrels: chr [1:4019] "S66366" "S48371" "S30713" "S43109" ...

pc_part2 <- pcairPartition(gds, kinobj = kingMat, kin.thresh = 2^(-3.5), div.thresh = -2^(-4.5), divobj = kingMat)
str(pc_part2)
# $ rels  : chr [1:955] "S48216" "S71706" "S25290" "S95760" ...
# $ unrels: chr [1:4016] "S66366" "S48371" "S30713" "S43109" ...

pc_part3 <- pcairPartition(gds, kinobj = kingMat, kin.thresh = 2^(-5.5), div.thresh = -2^(-4.5), divobj = kingMat)
str(pc_part3)
# $ rels  : chr [1:1008] "S26684" "S46646" "S89832" "S88925" ...
# $ unrels: chr [1:3963] "S66366" "S48371" "S30713" "S43109" ...


mypcair <- pcair(gds, kinobj = kingMat, kin.thresh = 2^(-4.5),
                 divobj = kingMat, snp.include = variant_id,
                 sample.include = sample_id, div.thresh = -2^(-4.5))
#Using kinobj and divobj to partition samples into unrelated and related sets
#Working with 4971 samples
#Identifying relatives for each sample using kinship threshold 0.0441941738241592
#Identifying pairs of divergent samples using divergence threshold -0.0441941738241592
#Partitioning samples into unrelated and related sets...
#Unrelated Set: 3994 Samples
#Related Set: 977 Samples
#Performing PCA on the Unrelated Set...
#Principal Component Analysis (PCA) on genotypes:
#Calculating allele counts/frequencies ...
#...



saveRDS(mypcair, paste0("6,18_1it", "pcair.rds"))


# plot 1st iteration PCs with script
Rscript pca_plots.R 6,18_1itpcair.rds --out_prefix 6,18_1it --phenotype_file annot.rds --group race_or_ethnicity

#Plot correaltion by PCs
qsub -q new.q calculate_snp_pc_corr.sh
qsub -q new.q plot_snp_pc_corr.sh #Not working yet... need to install Topmed Pipeline or SSH into Adrienne's folder



# 1st ititeration PC-Relate (restart R)
library(SeqArray)
library(GENESIS)
library(SeqVarTools)

gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("pruned_excludedRegions_andChr7.rds")
sample_id <- readRDS("keep_samples.rds")
mypcair <- readRDS("6,18_1itpcair.rds")

seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id)
## of selected samples: 4,971
## of selected variants: 171,679
seqData <- SeqVarData(gds)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
str(iterator)
#Formal class 'SeqVarBlockIterator' [package ""] with 6 slots
#  ..@ variantBlock : int 10000
#  ..@ variantFilter:List of 18
mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, seq(3)],
                    training.set = mypcair$unrels, sample.include = sample_id)
#4971 samples to be included in the analysis...
#Betas for 3 PC(s) will be calculated using 3994 samples in training.set...
#Running PC-Relate analysis for 4971 samples using 171679 SNPs in 18 blocks..
#...


saveRDS(mypcrel, paste0("6,18_1it", "pcr_obj.rds"))
pcr_mat <- pcrelateToMatrix(mypcrel, scaleKin = 1)
#Using 4971 samples provided
#Identifying clusters of relatives...
#    4971 relatives in 1 clusters; largest cluster = 4971
#Creating block matrices for clusters...
#0 samples with no relatives included
saveRDS(pcr_mat, paste0("6,18_1it", "pcr_mat.rds"))



# Plot 1st ititeration PC-Relate (restart R)
library(ggplot2)
library(SNPRelate)

rel <- readRDS("6,18_1itpcr_obj.rds")

kinship.df <- rel$kinBtwn

png(paste0("6,18", "1it_kinship.png"))
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
kingMat <- readRDS("6,18king_grm.rds")
pcrelate_matrix <- readRDS("6,18_1itpcr_mat.rds")

pc_part1 <- pcairPartition(gds, kinobj = pcrelate_matrix, kin.thresh = 2^(-4.5), div.thresh = -2^(-4.5), divobj = kingMat)
str(pc_part1)


mypcair <- pcair(gds, kinobj = pcrelate_matrix, kin.thresh = 2^(-4.5),
                 divobj = kingMat, snp.include = variant_id,
                 sample.include = sample_id, div.thresh = -2^(-4.5))

saveRDS(mypcair, paste0("6,18", "pcair.rds"))



# plot 1st iteration PCs with script
Rscript pca_plots.R 6,18pcair.rds --out_prefix 6,18 --phenotype_file annot.rds --group race_or_ethnicity

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
mypcair <- readRDS("6,18_1itpcair.rds")

seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id)
## of selected samples: 4,971
## of selected variants: 171,679
seqData <- SeqVarData(gds)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
#Formal class 'SeqVarBlockIterator' [package ""] with 6 slots
#  ..@ variantBlock : int 10000
#  ..@ variantFilter:List of 18

mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, seq(3)],
                    training.set = mypcair$unrels, sample.include = sample_id)
#4971 samples to be included in the analysis...
#Betas for 3 PC(s) will be calculated using 3994 samples in training.set...
#Running PC-Relate analysis for 4971 samples using 171679 SNPs in 18 blocks...
#...

saveRDS(mypcrel, paste0("6,18", "pcr_obj.rds"))
pcr_mat <- pcrelateToMatrix(mypcrel, thresh = 2^(-4.5), scaleKin = 2)
#Using 4971 samples provided
#Identifying clusters of relatives...
#    1892 relatives in 893 clusters; largest cluster = 26
#Creating block matrices for clusters...
#3079 samples with no relatives included
#Putting all samples together into one block diagonal matrix
saveRDS(pcr_mat, paste0("6,18", "pcr_mat.rds"))



# Plot 2nd ititeration PC-Relate
library(ggplot2)
library(SNPRelate)

rel <- readRDS("6,18pcr_obj.rds")

kinship.df <- rel$kinBtwn

png(paste0("6,18", "kinship.png"))
ggplot(kinship.df, aes_string("k0", "kin")) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) + #Note: if you get a "partial transparancy is not supported..." error remove "alpha" argument
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()

