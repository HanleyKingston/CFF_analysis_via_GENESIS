library(SeqArray)
library(SNPRelate)
library(argparser)
sessionInfo()

# read arguments
argp <- arg_parser("LD pruning")
argp <- add_argument(argp, "gds_file", help="GDS file")
argp <- add_argument(argp, "--out_file", help="output file name", default="pruned_snps.rds")
argp <- add_argument(argp, "--sample_id", help="RDS file with vector of sample.id to include")
argp <- add_argument(argp, "--variant_id", help="RDS file with vector of variant.id to include")
argp <- add_argument(argp, "--maf", help="minimum MAF for variants to include", default=0.05)
argp <- add_argument(argp, "--missing", help="maximum missing call rate for variants to include", default=0.05)
argp <- add_argument(argp, "--r_threshold", help="r threshold for LD", default=sqrt(0.1))
argp <- add_argument(argp, "--window_size", help="window size in Mb", default=10)
argp <- add_argument(argp, "--autosome_only", help = "whether to include autosomes only or also X-chromosome", default=TRUE)
argp <- add_argument(argp, "--build", help = "map build to filter problematic regions")
argp <- add_argument(argp, "--exclude_regions", help = "a dataframe of regions to exclude from pruned_SNPs, with colnames: chrom, start.base, end.base, comment ... note: only used if build is also provided")
argv <- parse_args(argp)
print(argv)

# parse file paths
gds.file <- argv$gds_file
out.file <- argv$out_file
sample.id <- if (!is.na(argv$sample_id)) readRDS(argv$sample_id) else NULL
variant.id <- if (!is.na(argv$variant_id)) readRDS(argv$variant_id) else NULL
exclude_regions <- if( !is.na(argv$exclude_regions)) readRDS(argv$exclude_regions) else NULL


# open GDS file
gds <- seqOpen(gds.file)
  
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
  
PCAcorr_snps <- seqGetData(gds, "variant.id")
saveRDS(PCAcorr_snps, "PCAcorr_snps.rds")


# run LD pruning
snpset <- snpgdsLDpruning(gds,
                          sample.id = sample.id,
                          snp.id = variant.id,
                          maf = argv$maf,
                          missing.rate = argv$missing,
                          method = "corr",
                          slide.max.bp = argv$window_size * 1e6, 
                          ld.threshold = argv$r_threshold,
                          autosome.only = argv$autosome_only
                          )

# convert list with one element per chrom to vector
pruned <- unlist(snpset, use.names=FALSE)


pruned <- intersect(pruned, PCAcorr_snps)

print(length(pruned))


saveRDS(pruned, file=out.file)
