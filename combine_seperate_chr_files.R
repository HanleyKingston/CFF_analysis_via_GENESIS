args <- commandArgs(trailingOnly = TRUE)

chr <- readRDS(paste0(args[1], "chr_1pruned.rds"))

for(i in 2:22){
chr_temp <- readRDS(paste0(args[1], "chr_", i, "pruned.rds"))
chr <- append(chr, chr_temp)
}

str(chr)

saveRDS(chr, paste0(args[1], "prunedSNPs.rds"))
