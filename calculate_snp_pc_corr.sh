#R -q --vanilla --args --pca-file 6,15pcair.rds --gds-file CFF_sid_onlyGT.gds --block-size 32768 --chromosome $SGE_TASK_ID --outfile snp_corr_chr${SGE_TASK_ID}.rds < calculate_snp_pc_corr.R

library(argparser)
library(tidyverse)
library(glue)
library(SeqVarTools)
library(SNPRelate)
library(dplyr)
library(tibble)

argp <- arg_parser("Correlation of variants with PCs")
argp <- add_argument(argp, "--pca-file", help="path to PCA file on disk")
argp <- add_argument(argp, "--gds-file", help="path to GDS file on disk")
argp <- add_argument(argp, "--chromosome", help="chromosome to compute correlations for")
argp <- add_argument(argp, "--outfile", help = "output file")
argp <- add_argument(argp, "--npcs", default = 3, help = "number of PCs to compute correlations for")
argp <- add_argument(argp, "--block-size", default = 10000, help = "genotype block size for computing correlations")
argv <- parse_args(argp)
print(argv)

gds_file <- argv$gds_file
pca_file <- argv$pca_file
n_pcs = argv$npcs
chromosome = argv$chromosome
block_size = argv$block_size
outfile = argv$outfile

# Read in PCair and subset to unrelated only.
pcair <- readRDS(pca_file)
pcs_all <- pcair$vectors
pcs_unrel <- pcs_all[rownames(pcs_all) %in% pcair$unrels, ]
dim(pcs_unrel)

# Map sample ids to ids in GDS Files.
key <- readr::read_tsv("/labdata12/CF_WGS2/cff_gwas/tables/key_cffwgs.tsv")
j <- match(rownames(pcs_unrel), key$sid)
rownames(pcs_unrel) <- key$vcf_id[j]


gds <- seqOpen(gds_file)

seqSetFilter(gds, sample.id = rownames(pcs_unrel))

# Additional variant filtering?
# Chr 1 for now
idx_chr <- seqGetData(gds, "chromosome") == chromosome
seqSetFilter(gds, variant.sel = idx_chr)

# Set up an iterator to iterate over all variants.
iterator <- SeqVarBlockIterator(gds, variantBlock=block_size)
n.iter <- length(variantFilter(iterator))
i <- 1
iterate <- TRUE
corr_list <- list()
while(iterate) {

  # Calculate correlations between genotype and PCs for this block.
  geno <- refDosage(iterator)
  variants <- data.frame(
    variant.id = seqGetData(iterator, "variant.id"),
    position = seqGetData(iterator, "position"),
    stringsAsFactors = FALSE
  )

  # Quick check, add code to fix if it fails.
  stopifnot(all.equal(rownames(geno), rownames(pcs_unrel)))

  corr <- cor(geno, pcs_unrel[, 1:n_pcs], use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    setNames(sprintf("PC", 1:n_pcs)) %>%
    rownames_to_column(var = "variant.id") %>%
    mutate(variant.id = as.integer(variant.id))
  corr_list[[i]] <- variants %>% left_join(corr, by = "variant.id")

  message(sprintf("Finished iteration %s/%s", i, n.iter))

  # Iterate.
  i <- i + 1
  iterate <- iterateFilter(iterator)

}
seqClose(gds)

corr <- bind_rows(corr_list) %>%
  mutate(chromosome = !!chromosome) %>%
  select(variant.id, chromosome, everything())

saveRDS(corr, file = outfile)
