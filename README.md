# CFF_analysis_via_GENESIS


## merge_ind_chr_files.R
Merge gds files (excluding X chromosome)

## Convert gds_ids (currently match VCF_ids from participants file) to pids from participants file
I still need to do this

## Save gds sample ids as a vector:
library(SeqArray)
gdsfile <- "/labdata12/CF_WGS2/Variants/CFF_5134_GDSs/seqArray_onlyGT/CFF_5134_chr21_onlyGT.gds"
gds <- seqOpen(gdsfile)
gds.id <- seqGetData(gds, "sample.id")
write(gds.id, "/home/hkings/DATA/sample_id_gds.txt")

## flag_exclude_in_participants.R
Trims participants pheotype file to excldue anyone not in gds, add a column for flagging or excluding indiviuals not in corresponding flag and exclude lists (also exclude duplicated samples)
(note: flag and exclude lsits were extracted manually from the info in sample_flag.txt

## Generate_sample_filter.R
This generates a vector of sample IDs to be used as a filter in LD pruning, PC and GRM generation, and association testing 

## Generate_variant_filter.R
This generates a list of variant IDs based on filter criteria... I generated a more and less stringent filter
both filters follows suggested filters from: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112 (see:[C] Hard-filter SNPs on multiple expressions using VariantFiltration)
my gds file already excldues the X chromosome

moderate filter also excludes:
 MAF < 0.01
 SNPs and indels that don't pass VQSR filters (index.snvPASS = FALSE | index.indelPASS = FALSE)

stringent filter also excludes:
  MAF < 0.05
 missingness by variant > 0.05
 indels
 anything not biallelic
 SNPs that don't pass VQSR filters (index.snvPASS = FALSE)

Note: can also filter by MAF and missingness in GENESIS's LD-pruning, but I chose to do it here so I can use the same filter for the association testing

## LD_prune.R



