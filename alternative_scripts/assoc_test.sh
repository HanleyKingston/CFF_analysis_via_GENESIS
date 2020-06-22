#This isn't working... not sure why. There is a script that does the same thing in the main repository

#$ -S /bin/sh
#$ -N assoc_test
#$ -cwd
#$ -j y
#$ -t 20-22

module load conda

R -q --vanilla --args CFF_sid_onlyGT.gds annot.rds 6_18pcr_mat.rds F508_count gaussian --out_file assoc_test${SGE_TASK_ID} --covars "PC1 PC2 PC3" --variant_id var_filter_SNVs_MAF0.05.rds --sample_id keep_samples.rds --chromosome $SGE_TASK_ID < assoc_test.R &

#Run with qsub -q new.q assoc_test.sh
