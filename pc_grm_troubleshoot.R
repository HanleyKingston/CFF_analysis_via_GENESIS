#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("Generate PCs and GRM") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("--out_prefix", help = "Prefix for output files",
               default = "") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs") %>%
  add_argument("--kin_thresh1", default = 5.5,
               help = "Kinship threshold for pcair (2 ^ -kin_thresh)") %>%
  add_argument("--div_thresh1", default = 5.5,
               help = "threshold for deciding if pairs are ancestrally divergent(-2 ^ -div_thresh)") %>%
  add_argument("--kin_thresh", default = 5.5,
               help = "Kinship threshold for pcair (2 ^ -kin_thresh)") %>%
  add_argument("--div_thresh", default = 5.5,
               help = "threshold for deciding if pairs are ancestrally divergent(-2 ^ -div_thresh)") %>%
  add_argument("--n_pcs", default = 3,
               "Number of PCs to pass to pcrelate") %>%
  add_argument("--keep_king", flag = TRUE, help = "Save KING-robust GRM")
argv <- parse_args(argp)

library(SeqArray)
library(GENESIS)
library(SeqVarTools)
library(SNPRelate)

sessionInfo()
print(argv)

if (!is.na(argv$variant_id)) {
  variant_id <- readRDS(argv$variant_id)
} else {
  variant_id <- NULL
}
if (!is.na(argv$variant_id)) {
  sample_id <- readRDS(argv$sample_id)
} else {
  sample_id <- NULL
}
kin_thresh1 <- 2 ^ (-argv$kin_thresh1)
div_thresh1 <- -2 ^ (-argv$div_thresh1)
kin_thresh <- 2 ^ (-argv$kin_thresh)
div_thresh <- -2 ^ (-argv$div_thresh)
out_prefix <- argv$out_prefix
gds <- seqOpen(argv$gds_file)

king <- snpgdsIBDKING(gds, snp.id = variant_id, sample.id = sample_id, type = "KING-robust")
kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id

if (argv$keep_king) {
  #Save king object
  saveRDS(king, paste0(out_prefix, "king_obj.rds"))
  #Save king matrix
  kingMat_temp <- kingMat * 2 # Scaled to match pc-relate GRM
  # coerces low values in matrix to 0
  kingMat_temp[kingMat_temp <= kin_thresh] <- 0
  saveRDS(kingMat_temp, paste0(out_prefix, "king_robust_grm.rds"))
  rm(kingMat_temp)
}


mypcair <- pcair(gds, kinobj = kingMat, kin.thresh = kin_thresh1, div.thresh = div_thresh1,
                 divobj = kingMat, snp.include = variant_id,
                 sample.include = sample_id)
print(str(mypcair))

#Save 1st iteration PCA object:
saveRDS(mypcair, paste0(out_prefix, "pcair_1it.rds"))

#Generate 1st iteration PC-Relate
seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id)
seqData <- SeqVarData(gds)
print("1st iteration PC-relate")
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, seq(argv$n_pcs)],
                    training.set = mypcair$unrels)
pcrelate_matrix <- pcrelateToMatrix(mypcrel, scaleKin=1)

#Save 1st iteration PC-Relate object:
saveRDS(mypcrel, paste0(out_prefix, "pcr_obj_1it.rds"))

#Generate 2nd iteration PC-Air:
pca <- pcair(gds, kinobj = pcrelate_matrix, kin.thresh = kin_thresh, div.thresh = div_thresh,
             divobj = kingMat, snp.include = variant_id,
             sample.include = sample_id)
print(str(pca))

saveRDS(pca, paste0(out_prefix, "pcair.rds"))

#Generate 2nd iteration PC-Relate object:
resetIterator(iterator, verbose = TRUE)
#iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)

print("2nd iteration PC-relate")
pcrel2 <- pcrelate(iterator, pcs = pca$vectors[, seq(argv$n_pcs)],
                   training.set = pca$unrels)

pcrelate_matrix2 <- pcrelateToMatrix(pcrel2, scaleKin = 2, thresh = kin_thresh)
saveRDS(pcrelate_matrix2, paste0(out_prefix, "pcr_grm.rds"))
saveRDS(pcrel2, paste0(out_prefix, "pcr_obj.rds"))

