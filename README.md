# CFF_analysis_via_GENESIS

## Contents:
### tables:
dd_key_expanded_cffwgs.tsv - Data dictionary for Key table (including all aliases) dd_samples_cffwgs.tsv - Data dictionary for Samples table key_expanded_cffwgs.tsv - Key table, including all aliases (Warning: contains eDWID) samples_cffwgs.tsv - Samples table dd_key_cffwgs.tsv - Data dictionary for Key table (simple version) dd_participants_cffwgs.tsv - Data dictionary for Participants table key_cffwgs.tsv - Key table, simplified version
participants_cffwgs.tsv - Participants table
tables_cffwgs.RData - R workspace containing all tables (with object types set)
### scripts


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
Generate a list of pruned SNPs to include in PC and GRM analyses

## PC_and_grm_script2.R
Arguments:
1. gds_file: the file path to the gds file (with .vcf.gds extension)
2. LD-pruning R object (a list of variants to incldue)
3. keep_variants: a file path to a list of variants to keep (must match corresponding rownames in phenotype and gds - saved as an R object
4. keep_samples: a file path to a list of samples to keep (must match corresponding smaple IDs in phenotype and gds files - default is familyID_SUBJID) - saved as an R object
5. text to uniquely identify plots and figures

ex. run with command: Rscript PC_and_grm_script2.R CFF_5134_onlyGT.gds "pruned.rds" "keep_var_stringent.rds" "keep_samples.rds" "LDsqrt0.1" #include "& > LDsqrt0.1_PCs_grm_script.out" to run concurrently with other processes and save output to a file (saving output only saves some basic info, I'm working on making it so it prints the whole console to file)




