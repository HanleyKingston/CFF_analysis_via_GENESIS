# LD-pruning

library(SeqArray)
library(SNPRelate)

# open GDS file
gds <- seqOpen("CFF_sid_onlyGT.gds") #This is the merged gds file ("CFF_5134_onlyGT.gds") from the seqArray_onlyGT folder, but I changed the sample.id's from the VCF_ids to the sids
variant_id <- readRDS("keep_var_stringent.rds")
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

# convert list with one element per chrom to vector
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)
#[1] 394140

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
## of selected variants: 111,797,708
}
  
PCAcorr_snps <- seqGetData(gds, "variant.id")
saveRDS(PCAcorr_snps, "PCAcorr_snps.rds")

pruned_excludedRegions_andChr7 <- intersect(pruned, PCAcorr_snps)
length(pruned_excludedRegions_andChr7)
#[1] 380171

saveRDS(pruned_excludedRegions_andChr7, "pruned_excludedRegions_andChr7.rds")




# King
gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("pruned_excludedRegions_andChr7.rds")
sample_id <- readRDS("keep_samples.rds")

king <- snpgdsIBDKING(gds, snp.id = variant_id, sample.id = sample_id,
                      type = "KING-robust")

kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id

saveRDS(king, paste0("6,15", "king_obj.rds"))
saveRDS(kingMat, paste0("6,15", "king_grm.rds")) #Note: should save a sparse matrix if using for association testing but not to give to PC-AiR



# Plot King
library(ggplot2)
library(SNPRelate)

rel <- readRDS("6,15king_obj.rds")

kinship.df <- snpgdsIBDSelection(rel)

png(paste0("6,15", "_kinship.png"))
ggplot(kinship.df, aes_string("IBS0", "kinship")) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) + #Note: if you get a "partial transparancy is not supported..." error remove "alpha" argument
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()



# 1st itertion PC-Air
library(SeqArray)
library(GENESIS)

gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("pruned_excludedRegions_andChr7.rds")
sample_id <- readRDS("keep_samples.rds")
kingMat <- readRDS("6,15king_grm.rds")

mypcair <- pcair(gds, kinobj = kingMat, kin.thresh = 2^(-3),
                 divobj = kingMat, snp.include = variant_id,
                 sample.include = sample_id, div.thresh = -2^(-4.5))

saveRDS(mypcair, paste0("6,15", "pcair_1it.rds"))



# 1st ititeration PC-Relate
gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("pruned_excludedRegions_andChr7.rds")
sample_id <- readRDS("keep_samples.rds")
mypcair <- readRDS("6,15pcair_1it.rds")

seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id)
seqData <- SeqVarData(gds)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, seq(12)],
                    training.set = mypcair$unrels, sample.include = sample_id)

saveRDS(mypcrel, paste0("6,15", "pcr_obj_1it.rds"))
pcr_mat <- pcrelateToMatrix(mypcrel, scaleKin = 1)
saveRDS(pcr_mat, paste0(argv$out_prefix, "pcr_mat.rds"))




# 1st ititeration PC-Relate
library(ggplot2)
library(SNPRelate)

rel <- readRDS("6,15pcr_obj_1it.rds")

kinship.df <- rel$kinBtwn

png(paste0("6,15", "1it_kinship.png"))
ggplot(kinship.df, aes_string("k0", "kin")) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color = "grey") +
    geom_point(alpha=0.2) + #Note: if you get a "partial transparancy is not supported..." error remove "alpha" argument
    ylab("kinship estimate") +
    ggtitle("kinship")
dev.off()



# 2nd itertion PC-Air
library(SeqArray)
library(GENESIS)

gds <- seqOpen("CFF_sid_onlyGT.gds")
variant_id <- readRDS("pruned_excludedRegions_andChr7.rds")
sample_id <- readRDS("keep_samples.rds")
kingMat <- readRDS("6,15king_grm.rds")
kingMat <- readRDS("6,15pcr_obj_1it.rds")


mypcair <- pcair(gds, kinobj = kingMat, kin.thresh = 2^(-3),
                 divobj = kingMat, snp.include = variant_id,
                 sample.include = sample_id, div.thresh = -2^(-4.5))

saveRDS(mypcair, paste0("6,15", "pcair_1it.rds"))



