# CFF_analysis_via_GENESIS

## Contents:
### tables:
samples_flag.txt - samples to flag (extracted from samples_flag_or_exclude.txt)
samples_excldue.txt - sampels to exclude (extracted from samples_flag_or_exclude.txt)
participants_cffwgs.tsv - participant phenotype data (note: not uploaded yet, need to make directory private first)
samples_flag_or_exclude.txt - a list of samples to flag or exclude (note:not uploaded yet)
keep_samples.rds - a vector of sample IDs to include in PC/GRM creation and assoc. testing
keep_var_stringent.rds - vector of variant IDs to include in LD-pruning, PC/GRM creation, and association testing... based on stringent QC and pruning criteria (see below)
  
### scripts
flag_exclude_in_particiapnts.R - annotate participants phenotype data based on sampels to flag or excldue
Generate_sample_filter.R - script to create a vector of sample IDs to include in PC/GRM creation and assoc. testing
Generate_variant_filter.R - script to create a vector of variant IDs to include in LD-pruning, PC/GRM creation, and association testing
LD_prune.R - generate a vector of SNPs pruned by LD to include in PC and GRM creation
PC_and_grm_script2.R - generate PCs and GRM through 2 iteratons of PCair and PCrelate (and plot first 3 PCs and kinship)



## merge_ind_chr_files.R
Merge gds files (excluding X chromosome)

## Save gds sample ids as a vector:
### in R:
library(SeqArray)  
gds <- seqOpen("CFF_5134_onlyGT.gds")  
gds.id <- seqGetData(gds, "sample.id")  
saveRDS(gds.id, "sample_id_gds.rds")  
seqClose(gds)

## flag_exclude_in_participants.R
Trims participants pheotype file to excldue anyone not in gds, add a column for flagging or excluding indiviuals not in corresponding flag and exclude lists (also exclude duplicated samples)
(note: flag and exclude lsits were extracted manually from the info in sample_flag.txt

## change_gds_ids.R
Convert gds_ids (currently match VCF_ids from participants file) to pids from participants file

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
Takes Arguments:
1. gds_file: the file path to the gds file (with .vcf.gds extension)
2. threshold for LD-pruning (given as the correlation value which should be the square route of R^2) - variants above the threshold (ie. in greater LD, are pruned)
3. keep_variants: a file path to a list of variants to keep (must match corresponding rownames in phenotype and gds - saved as an R object
4. keep_samples: a file path to a list of samples to keep (must match corresponding smaple IDs in phenotype and gds files - default is familyID_SUBJID) - saved as an R object

#ex. run with command: Rscript LD_prune.R CFF_5134_onlyGT.gds 0.316227766 "keep_variants.rds" "keep_samples.rds" & > LDsqrt0.1_PCs_grm_script.out

## PC_and_grm_script2.R
Arguments:
1. gds_file: the file path to the gds file (with .vcf.gds extension)
2. LD-pruning R object (a list of variants to incldue)
3. keep_variants: a file path to a list of variants to keep (must match corresponding rownames in phenotype and gds - saved as an R object
4. keep_samples: a file path to a list of samples to keep (must match corresponding smaple IDs in phenotype and gds files - default is familyID_SUBJID) - saved as an R object
5. text to uniquely identify plots and figures

ex. run with command: Rscript PC_and_grm_script2.R CFF_5134_onlyGT.gds "pruned.rds" "keep_var_stringent.rds" "keep_samples.rds" "LDsqrt0.1" #include "& > LDsqrt0.1_PCs_grm_script.out" to run concurrently with other processes and save output to a file (saving output only saves some basic info, I'm working on making it so it prints the whole console to file)

## PC_and_GRM_plots.R
Generate PCs plots, percent variance explained (scree) plots, and relatedness plots



