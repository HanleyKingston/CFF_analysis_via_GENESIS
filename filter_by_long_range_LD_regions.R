#' @param build Genome build to use when identifying regions to exclude from PCA because of high correlation (HLA, LCT, inversions)
#' @rdname filterVariants
#'
#' @import SeqArray
#' @importFrom utils data
#' @export
filterByPCAcorr <- function(gds, build="hg19", verbose=TRUE) {
    filt <- get(data(list=paste("pcaSnpFilters", build, sep="."), package="GWASTools"))
    chrom <- seqGetData(gds, "chromosome")
    pos <- seqGetData(gds, "position")
    pca.filt <- rep(TRUE, length(chrom))
    for (f in 1:nrow(filt)) {
        pca.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos & pos < filt$end.base[f]] <- FALSE
    }
    seqSetFilter(gds, variant.sel=pca.filt, action="intersect", verbose=verbose)
}
