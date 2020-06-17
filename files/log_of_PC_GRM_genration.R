# LD-pruning

library(SeqArray)
library(SNPRelate)
library(argparser)
sessionInfo()

# open GDS file
gds <- seqOpen(gds.file)

# run LD pruning
snpset <- snpgdsLDpruning(gds,
                          sample.id = keep_samples.rds,
                          snp.id = keep_var_stringent.rds,
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
