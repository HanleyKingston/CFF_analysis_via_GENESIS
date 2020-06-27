library(SeqArray)
library(SNPRelate)
library(argparser)
sessionInfo()

# read arguments
argp <- arg_parser("LD pruning")
argp <- add_argument(argp, "gds_file", help="GDS file")
argp <- add_argument(argp, "--out_prefix", help="output file name", default="")
argp <- add_argument(argp, "--sample_id", help="RDS file with vector of sample.id to include")
argp <- add_argument(argp, "--variant_id", help="RDS file with vector of variant.id to include")
argp <- add_argument(argp, "--maf", help="minimum MAF for variants to include", default=0.05)
argp <- add_argument(argp, "--missing", help="maximum missing call rate for variants to include", default=0.05)
argp <- add_argument(argp, "--r_threshold", help="r threshold for LD", default=sqrt(0.1))
argp <- add_argument(argp, "--window_size", help="window size in Mb", default=10)
argp <- add_argument(argp, "--autosome_only", help = "whether to include autosomes only or also X-chromosome", default=TRUE)
argp < -add_argument(argp, "--chromosome", help = "chromosome number")
argv <- parse_args(argp)
print(argv)

# parse file paths
gds.file <- argv$gds_file
out.file <- argv$out_file
sample.id <- if (!is.na(argv$sample_id)) readRDS(argv$sample_id) else NULL
variant.id <- if (!is.na(argv$variant_id)) readRDS(argv$variant_id) else NULL


# open GDS file
gds <- seqOpen(gds.file)
  
if (!ia.na(argv$chromosome) {
  seqSetFilterChrom(gds, chromosome)
  chrom_vars <- seqGetData(gds, "variant.id")
  variant.id <- intersect(variant.id, chrom_vars)
  seqResetFilter(gds)
}
    
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

print(length(pruned))

saveRDS(pruned, paste0(argv$out_prefix, "pruned.rds"))
