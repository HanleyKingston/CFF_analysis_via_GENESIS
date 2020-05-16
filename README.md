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


