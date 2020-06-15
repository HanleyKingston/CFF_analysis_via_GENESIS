#CFF_sid_onlyGT.gds

library(SeqArray)
library(utils)

argp <- arg_parser("Run KING-robust") %>%
  add_argument("gds_file", help = "GDS file") %>%


gds <- argv$gds_file

filterByPCAcorr <- function(gds, build="hg38", verbose=TRUE) {
    filt <- get(data(list=paste("pcaSnpFilters", build, sep="."), package="GWASTools"))
    chrom <- seqGetData(gds, "chromosome")
    pos <- seqGetData(gds, "position")
    pca.filt <- rep(TRUE, length(chrom))
    for (f in 1:nrow(filt)) {
        pca.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos & pos < filt$end.base[f]] <- FALSE
    }
    seqSetFilter(gds, variant.sel=pca.filt, action="intersect", verbose=TRUE)
}
