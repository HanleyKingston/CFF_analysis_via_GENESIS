add_regions <- data.frame(t(c(7, 1, 159345973, "exclude_chr7_to_exclude_deltaF508")))
colnames(add_regions) <- colnames(filt)
rownames(add_regions) <- "chr7"

print(as.data.frame(add_regions))

saveRDS(as.data.frame(add_regions), "exclude_regions_chr7.rds")
