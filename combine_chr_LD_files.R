args <- commandArgs(trailingOnly = TRUE)

chr <- readRDS(paste0(args[1], "chr1_pruned.rds"))

for(i in 2:22){
chr_temp <- readRDS(paste0(args[1], "chr_", i, "pruned.rds"))
chr <- append(chr, chr_temp)
}

 str(chr)
 #int [1:171975] 1683 2190 2215 3755 6388 10751 10852 11502 13043 13086 ...

saveRDS(chr, paste0(args[1], "6_26_prunedSNPs.rds"))
