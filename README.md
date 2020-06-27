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

## change_gds_ids.R
1. Convert gds_ids (currently match VCF_ids from participants file) to sids from based on sample_names_key.txt (taken from SampleDropandFlag.tsv)
2. Save gds IDs as a vector

## generate_sample_filter_and_phenotype_df.R
1. Trims participants pheotype file to excldue anyone not in gds
2. Add columns:
include_in_analysis: marks indiviuals "exclude" based on SampleDropandFlag file.
(note: flag and exclude lsits were extracted manually from the info in sample_flag.txt
site: based on VCF_ID (likely to have some errors due to duplicate samples)
F508_count: based on CFTR genotype column
race_or_ethnicity: extracted based on binary race and ethnicity measures (note: hispanic white people are labeled hispanic and Native American white people are labeled Native American. Anyone else with mutliple races or ethnicities is labeled admixed_or_other

## Generate_varaint_filters_from_SNPs_filtered.R
Check that varaint filter variant IDs match teh GDS variant IDs by chromosome and position
and generate variant filters for LD-pruning and association-testing steps

## exclude_regions_beforeLD_prune.R
This will create a dataframe of all of chr7 to exclude from pruned SNP list for PC and GRM generation

## ld_pruning.sh
### sh ld_pruning.sh
runs ld_pruning.R with the following arguments:
R -q --vanilla --args CFF_sid_onlyGT.gds --sample_id keep_samples.rds --variant_id pre_LD_SNP_filter.rds --maf 0.05 --missing 0.05 --window_size 1 --autosome_only TRUE < ld_pruning.R > 6_26ld_pruning.log &

generate a vector of SNPs pruned by LD to include in PC and GRM creation

Takes Arguments:
1. req: gds_file: .gds file (with a chracter vector of sid as sample IDs)
2. opt: out_file (default="pruned_snps.rds")
3. opt: sample_id: a character vector of samples to keep (must match corresponding smaple IDs in phenotype and gds files) - saved as an R object
4. opt: variant_id: a vector of variants to keep (must match corresponding rownames in phenotype and gds) - saved as an R object
5. opt: maf: minimum MAF for variants to include (default=0.05)
6. opt: missing: maximum missing call rate for variants to include, (default=0.05)
7. opt: threshold: threshold for LD-pruning (given as the correlation value which should be the square route of R^2) - variants above the threshold (ie. in greater LD, are pruned) (default=sqrt(0.1))
8. opt: window_size (default = 1) (slide.max.bp = argv$window_size * 1e6)

### Rscript combine_chr_LD_files.R


## king_grm.R
### R -q --vanilla --args CFF_sid_onlyGT.gds --out_prefix 6_26 --variant_id 6_26_prunedSNPs.rds --sample_id keep_samples.rds --autosome_only TRUE < king_grm.R > 6_26king_grm.log &

## kinship plots.R
### Rscript plot_kinship.R 6_26king_out.rds --is_king --out_prefix 6_26_king

## pcair.R
#To determine a kin_thresh:
library(GENESIS)
library(SeqArray
gds <- seqOpen("CFF_sid_onlyGT.gds")
kingMat <- readRDS("6_25king_grm.rds")
pc_part <- pcairPartition(gds, kinobj = kingMat, kin.thresh = 2^(-4.5), div.thresh = -2^(-4.5), divobj = kingMat)
str(pc_part)
### R -q --vanilla --args CFF_sid_onlyGT.gds 6_26king_grm.rds 6_26king_grm.rds --out_prefix 6_26_1it --variant_id 6_26_prunedSNPs.rds --sample_id keep_samples.rds --kin_thresh 0.044194 --div_thresh -0.044194 < pcair.R > 6_26_1itpc_air.log &
#0.044194 = 2^(-9/2)
#0.0625 = 2^(-4)
#0.125 = 2^(-3)

#Fast way to determine how many PCs to include (note: should also look at PC correlation plots):
pca <- readRDS("6_26_1itpcair.rds")
plot(seq(12),100*pca$varprop[1:12])


## pcrelate.R
### R -q --vanilla --args CFF_sid_onlyGT.gds 6_26_1itpcair.rds --out_prefix 6_26_1it --n_pcs 4 --variant_id 6_26_prunedSNPs.rds --sample_id keep_samples.rds --scale_kin 1 --small_samp_correct --variant_block 100000 < pcrelate.R > 6_26_1itpcrelate.log &

## kinship plots.R
### Rscript plot_kinship.R 6_26_1itpcrelate.rds --out_prefix 6_26_1it_PC-rel


## pcair.R
### R -q --vanilla --args CFF_sid_onlyGT.gds 6_26pcr_mat.rds 6_26king_grm.rds --out_prefix 6_26 --variant_id 6_26_prunedSNPs.rds --sample_id keep_samples.rds --kin_thresh 0.044194 --div_thresh -0.044194 < pcair.R > 6_26pc_air.log &
#0.044194 = 2^(-9/2)

## pcrelate.R
### R -q --vanilla --args CFF_sid_onlyGT.gds 6_26_1itpcair.rds --out_prefix 6_26 --n_pcs 4 --variant_id 6_26_prunedSNPs.rds --sample_id keep_samples.rds --scale_kin 1 --small_samp_correct --variant_block 100000 < pcrelate.R > 6_26pcrelate.log &

## kinship plots.R
### Rscript plot_kinship.R 6_26_pcrelate.rds --out_prefix 6_26_PC-rel

Arguments for PC-Relate and PC-Air (needs updating):
1. req: gds_file: .gds file (with a chracter vector of sid as sample IDs)
2. opt : out_prefix
3. opt: variant_id: a vector of variants to keep (must match corresponding rownames in phenotype and gds) - use pruned_snps.rds from LD pruning step - saved as an R object
4. opt: sample_id: a character vector of samples to keep (must match corresponding smaple IDs in phenotype and gds files) - saved as an R object2. LD-pruning R object (a list of variants to incldue)
5. opt: kin_thresh: Kinship threshold for pcair (2 ^ -kin_thresh) (default = 5.5)
6. opt: div_thresh: threshold for deciding if pairs are ancestrally divergent (-2 ^ -kin_thresh) (default = 5.5)
7. opt: n_pcs: number of PCs to pass to PC-Relate (default = 3)
8. opt: text to uniquely identify plots and figures
9. opt: keep_king: if passed TRUE, will also save the GRM from KING robust


## add_phenotype_identifiers_to_kinship_obj.R 
-Needs to be fixed to run efficiently-
To color kinship plot based on ancestry, study, etc.
Arguments:
1. phenotype data frame as .rds object
2. a character vector of samples to keep (as a .rds object)
3. kinship object (not matrix) from King (as a .rds object)
4. The type of relatedness object: either "King" or "PC_Relate"
### R -q --vanilla --args phenotype.rds keep_samples.rds higherLDking_obj.rds King < add_phenotype_identifiers_to_kinship_obj.R & > record.txt
*This takes hours to run, so include &

## Generate_annotated_phenotype_df.R
Add PCs to phenotype data and produce an annotated dataframe to be used in pca_plots.R and assoc_test.R

## pca_plots.R
Plots a scree plot (percent variance explained), cord plot, and pairwise PC comparisons to further anylize PCs. and phenotype file as an annotated data frame. Takes PC-AiR and PC-Relate .rds objects and phenotype file as an annotated data frame
### Rscript pca_plots.R 6_26pcair.rds --out_prefix 6_26pcair --phenotype_file annot.rds --group race_or_ethnicity


## Exclude_identical_twin_from_samples.R
Create a new sample filter that excludes idenitcal twins to be used in assoc_test.R

## assoc/assoc_test.sh
### sh assoc_test.sh
Runs assoc_test.R with the following arguments:
CFF_sid_onlyGT.gds annot.rds 6_18pcr_mat.rds F508_count gaussian --out_prefix "$f" --covars "PC1 PC2 PC3" --variant_id var_filter_SNVs_MAF0.05.rds --sample_id keep_samples.rds --chromosome "$f" < assoc_test2.R &
1. req: gds_file: .gds file (with a chracter vector of sid as sample IDs)
2. req: phenotype file (as an annotated dataframe R object)
3. req: outcome
4. req: faily (either binomial, gaussian (continuous), or poisson (discrete numeric))
5. opt: out_prefix (to parallelize, this should be the chrosomome name)
6. opt: covariates (as a string, each covariate should be seperated by a space)
7. opt: variant_id: a vector of variants to keep (must match corresponding rownames in phenotype and gds) - use pruned_snps.rds from LD pruning step - saved as an R object
5. opt: sample_id: a character vector of samples to keep (must match corresponding smaple IDs in phenotype and gds files) - saved as an R object2. LD-pruning R object (a list of variants to incldue)
7. opt: chromosome to fitler by

## combine_chr_assoc_files.R
recombines association test files into one file (assoc.rds) - chromosomes will be in numerical order
Takes as an argument:
1. prefix for input and output (input files should be in the form: "prefix"chr1assoc.rds)
### Rscript combine_chr_assoc_files.R F508del
### Rscript combine_chr_assoc_files.R sex

## assoc/assoc_plots.R
###  R -q --vanilla --args F508del_assoc.rds --out_prefix CFF_F508 < assoc_plots.R &
###  R -q --vanilla --args sex_assoc.rds --out_prefix sex < assoc_plots.R &

