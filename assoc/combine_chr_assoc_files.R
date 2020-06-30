args <- commandArgs(trailingOnly = TRUE)

chr <- readRDS(paste0(args[1], "_chr_1assoc.rds"))

for(i in 2:22){
chr_temp <- readRDS(paste0(args[1], "_chr_", i, "assoc.rds"))
chr <- rbind(chr, chr_temp)
}

saveRDS(chr, paste0(args[1], "assoc.rds"))
