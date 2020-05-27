# CFF_analysis_via_GENESIS

## Contents:
### files and tables (not all of these are uploaded but all are needed to run analysis)
samples_flag.txt - samples to flag (extracted from samples_flag_or_exclude.txt)
samples_excldue.txt - sampels to exclude (extracted from samples_flag_or_exclude.txt)
participants_cffwgs.tsv - participant phenotype data (note: not uploaded yet, need to make directory private first)
samples_flag_or_exclude.txt - a list of samples to flag or exclude (note:not uploaded yet)
CFF_sid_onlyGT.gds - gds file with sid as sample IDs
keep_samples.rds - a vector of sample IDs to include in PC/GRM creation and association testing
keep_var_stringent.rds - vector of variant IDs to include in LD-pruning, PC/GRM creation, and association testing... based on stringent QC and pruning criteria (see below)
pc_grm_output.txt - manually saved command line output from generating PCs and GRMs
keep_sample_noTwins.rds - a vector of sample IDs to include in association testing
  
### scripts
flag_exclude_in_particiapnts.R - annotate participants phenotype data based on sampels to flag or excldue
Generate_sample_filter.R - script to create a vector of sample IDs to include in PC/GRM creation and assoc. testing
Generate_variant_filter.R - script to create a vector of variant IDs to include in LD-pruning, PC/GRM creation, and association testing
ld_pruning.R - generate a vector of SNPs pruned by LD to include in PC and GRM creation
pcs_and_grm.R - generate PCs and GRM through 2 iteratons of PCair and PCrelate (and plot first 3 PCs and kinship)
PC_andGRM_plots.R - plots PCs w/ more features, simple plot of percent variance explained and relatedness plot. Takes PC-AiR and PC-Relate .rds objects
Generate_annotated_phenotype_df.R - Add PCs to phenotype data and produce an annotated dataframe to be used in pca_plots.R and assoc_test.R
pca_plots.R - plots a scree plot (percent variance explained), cord plot, and pairwise PC comparisons to further anylize PCs
Exclude_identical_twin_from_samples.R - Create a new sample filter that excludes idenitcal twins to be used in assoc_test.R
assoc_test.R


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

## ld_pruning.R
generate a vector of SNPs pruned by LD to include in PC and GRM creation

Takes Arguments:
1. req: gds_file: .gds file (with a chracter vector of sid as sample IDs)
2. opt: out_file (default="pruned_snps.rds")
3. opt: sample_id: a cahracter vector of samples to keep (must match corresponding smaple IDs in phenotype and gds files) - saved as an R object
4. opt: variant_id: a vector of variants to keep (must match corresponding rownames in phenotype and gds) - saved as an R object
5. opt: maf: minimum MAF for variants to include (default=0.05)
6. opt: missing: maximum missing call rate for variants to include, (default=0.05)
7. opt: threshold: threshold for LD-pruning (given as the correlation value which should be the square route of R^2) - variants above the threshold (ie. in greater LD, are pruned) (default=sqrt(0.1))
8. opt: window_size (default = 1)

### script ld_pruning.R CFF_sid_onlyGT.gds --sample_id keep_samples.rds --variant_id keep_var_stringent.rds --window_size 1


## pcs_and_grm.R
generate PCs and GRM through 2 iteratons of PCair and PCrelate (and plot first 3 PCs and kinship)

Arguments:
1. req: gds_file: .gds file (with a chracter vector of sid as sample IDs)
2. opt : out_prefix
3. opt: variant_id: a vector of variants to keep (must match corresponding rownames in phenotype and gds) - use pruned_snps.rds from LD pruning step - saved as an R object
4. opt: sample_id: a cahracter vector of samples to keep (must match corresponding smaple IDs in phenotype and gds files) - saved as an R object2. LD-pruning R object (a list of variants to incldue)
5. opt: kin_thresh: Kinship threshold for pcair (2 ^ -kin_thresh) (default = 5.5)
6. opt: div_thresh: threshold for deciding if pairs are ancestrally divergent (-2 ^ -kin_thresh) (default = 5.5)
7. opt: n_pcs: number of PCs to pass to PC-Relate (default = 3)
8. opt: text to uniquely identify plots and figures
9. opt: keep_king: if passed TRUE, will also save the GRM from KING robust

### Rscript pcs_and_grm.R CFF_sid_onlyGT.gds --out_prefix CFF_LDsqrt0.1 --variant_id pruned_snps.rds --sample_id keep_samples.rds --kin_thresh 4.5 --div_thresh 4.5
#include "& > LDsqrt0.1_PCs_grm_script.out" to run concurrently with other processes and save output to a file (saving output only saves some basic info, I'm working on making it so it prints the whole console to file)


## PC_and_GRM_plots.R
plots PCs w/ more features, simple plot of percent variance explained and relatedness plot.

## Generate_annotated_phenotype_df.R
Add PCs to phenotype data and produce an annotated dataframe to be used in pca_plots.R and assoc_test.R

## pca_plots.R
Plots a scree plot (percent variance explained), cord plot, and pairwise PC comparisons to further anylize PCs. and phenotype file as an annotated data frame. Takes PC-AiR and PC-Relate .rds objects and phenotype file as an annotated data frame

## Exclude_identical_twin_from_samples.R
Create a new sample filter that excludes idenitcal twins to be used in assoc_test.R
assoc_test.R

## assoc_test.R



