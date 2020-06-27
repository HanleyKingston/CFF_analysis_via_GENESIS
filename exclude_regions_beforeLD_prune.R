#Make a dataframe of chr 7 SNPs to incldue in filter of SNPs to exclude in LD-steps
add_regions <- data.frame(t(c(7, 1, 159345973, "exclude_chr7_to_exclude_deltaF508")))
colnames(exclude_regions) <- c("chrom", "start.base", "end.base", "comment")
rownames(exclude_regions) <- "chr7"
exclude_regions <- as.data.frame(exclude_regions)

saveRDS(exclude_regions, "exclude_regions_chr7.rds")

#Filter out problematic LD regions
if( !is.na(argv$build)){ #Need to check what this will do if not build is supplied
    filt <- get(data(list=paste("pcaSnpFilters", argv$build, sep="."), package="GWASTools"))
    filt <- rbind(filt, exclude_regions)
    chrom <- seqGetData(gds, "chromosome")
    pos <- seqGetData(gds, "position")
    pca.filt <- rep(TRUE, length(chrom))
    for (f in 1:nrow(filt)) {
        pca.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos & pos < filt$end.base[f]] <- FALSE
    }
    seqSetFilter(gds, variant.sel=pca.filt, action="intersect", verbose=TRUE)
}

non_prob_snps <- seqGetData(gds, "variant.id")
saveRDS(non_prob_snps, "non_prob_snps.rds")

#Save final filter to pass to LD-pruning
pre_LD_SNP_filter <- intersect(non_prob_snps, variant.id)
saveRDS(pre_LD_SNP_filter, "pre_LD_SNP_filter.rds")
