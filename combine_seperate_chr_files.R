chr <- readRDS("chr1assoc.rds")

for(i in 2:22){
chr_temp <- readRDS(paste0("chr", i, "assoc.rds"))
chr <- rbind(chr, chr_temp)
}

saveRDS(chr, "assoc.rds")
