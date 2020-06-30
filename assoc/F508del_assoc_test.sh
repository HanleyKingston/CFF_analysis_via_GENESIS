#$ -S /bin/sh
#$ -N F508del_assoc_test
#$ -cwd
#$ -j y
#$ -t 1-22

for f in {1..22}; do
echo "Assoc_testing '$f' "
R -q --vanilla --args CFF_sid_onlyGT.gds annot.rds 6_26pcr_mat.rds F508_count poisson --out_prefix F508del_chr_${SGE_TASK_ID} --covars "PC1 PC2 PC3 PC4 site" --variant_id SNPS_bi_GATK_VQSR_MAF0.05_miss0.05.rds --sample_id keep_samples_noTwins.rds --chromosome $SGE_TASK_ID < assoc_test.R > F508del_assoc_test${SGE_TASK_ID}.log &
sleep 30
done

#qsub -q new.q F508del_assoc_test.sh
