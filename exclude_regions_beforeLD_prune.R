add_regions <- data.frame(t(c(7, 1, 159,345,973, "to exclude deltaF508")))
colnames(add_regions) <- colnames(filt)
rownames(add_regions) <- "chr7"
saveRDS(add_regions, "exclude_regions_chr7.rds")
