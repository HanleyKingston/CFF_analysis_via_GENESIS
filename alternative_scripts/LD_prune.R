#Arguments:
#1. gds_file: the file path to the gds file (with .vcf.gds extension)
#2. threshold for LD-pruning (given as the correlation value which should be the square route of R^2) - variants above the threshold (ie. in greater LD, are pruned)
#3. keep_variants: a file path to a list of variants to keep (must match corresponding rownames in phenotype and gds - saved as an R object
#4. keep_samples: a file path to a list of samples to keep (must match corresponding smaple IDs in phenotype and gds files - default is familyID_SUBJID) - saved as an R object

#To test: args <- c("CFF_sid_onlyGT.gds", "0.316227766", "keep_var_stringent.rds", "keep_samples.rds")
#ex. run with command: Rscript LD_prune.R CFF_pid.gds 0.316227766 "keep_var_stringent.rds" "keep_samples.rds" &

args <- commandArgs(trailingOnly = TRUE)
print(args)

library(SeqArray)
library(SNPRelate)

##Open GDS file:
gdsfmt::showfile.gds(closeall=TRUE) #makes sure gds file isn't already open
gds <- seqOpen(args[1])
#Note: can access gds with: seqGetData(gds, "___")

thresh <- as.numeric(args[2])
keep_vars <- readRDS(file = args[3])
keep_samples <- readRDS(file = args[4])

class(seqGetData(gds, "sample.id"))
keep_samples <- as.vector(factor(keep_samples))

print(paste("LD applied at threshold:", thresh))
set.seed(100) # LD pruning has a random element; so make this reproducible
# LD pruning to get variant set:
snpset <- snpgdsLDpruning(gds, method="corr", ld.threshold=thresh, slide.max.bp = 1e6, num.thread = 25,  snp.id = keep_vars, sample.id = keep_samples)
