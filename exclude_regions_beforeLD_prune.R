library(SeqArray)
library(GENESIS)

#Make a dataframe of chr 7 SNPs to incldue in filter of SNPs to exclude in LD-steps
exclude_regions <- data.frame(t(c(7, 1, 159345973, "exclude_chr7_to_exclude_deltaF508")))
colnames(exclude_regions) <- c("chrom", "start.base", "end.base", "comment")
rownames(exclude_regions) <- "chr7"
exclude_regions <- as.data.frame(exclude_regions)

saveRDS(exclude_regions, "exclude_regions_chr7.rds")


#Filter out problematic LD regions
filt <- get(data(list=paste("pcaSnpFilters", "hg38", sep="."), package="GWASTools"))
filt <- rbind(filt, exclude_regions)
gds <- seqOpen("CFF_sid_onlyGT.gds")
chrom <- seqGetData(gds, "chromosome")
pos <- seqGetData(gds, "position")
pca.filt <- rep(TRUE, length(chrom))
for (f in 1:nrow(filt)) {
    pca.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos & pos < filt$end.base[f]] <- FALSE
}
seqSetFilter(gds, variant.sel=pca.filt, action="intersect", verbose=TRUE)


non_prob_snps <- seqGetData(gds, "variant.id")
saveRDS(non_prob_snps, "non_prob_snps.rds")

#Read in variant ID's that passed QC:
variant.id <- readRDS("SNPS_bi_GATK_VQSR.rds")
length(variant.id)
#[1] 84360856


#Save final filter to pass to LD-pruning
pre_LD_SNP_filter <- intersect(non_prob_snps, variant.id)
length(pre_LD_SNP_filter)
#[1] 81569276
saveRDS(pre_LD_SNP_filter, "pre_LD_SNP_filter.rds")
