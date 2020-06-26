#$ -S /bin/sh
#$ -N calculate_snp_pc_corr
#$ -cwd
#$ -j y
#$ -t 1-22

module load conda

R -q --no-save --args --pca-file /pc_correlations6,15/6,15pcair.rds --gds-file CFF_sid_onlyGT.gds --block-size 32768 --chromosome $SGE_TASK_ID --outfile /pc_correlations6,15/snp_corr_chr${SGE_TASK_ID}.rds < calculate_snp_pc_corr.R



#Run with qsub -q new.q -v R_LIBS=/home/amstilp/devel/analysis_pipeline_cff_wgs/R_library calculate_snp_pc_corr.sh
