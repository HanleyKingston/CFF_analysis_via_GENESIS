for f in {1..22}; do
echo "Assoc_testing '$f' "
R -q --vanilla --args CFF_sid_onlyGT.gds annot.rds 6-18pcr_mat.rds F508_count gaussian --out_prefix "$f" --covars "PC1 PC2 PC3" --variant_id var_filter_SNVs_MAF0.05.rds --sample_id keep_samples.rds --chromosome "$f" < assoc_test.R &
sleep 30
done

#assoc_test.sh


#Alt. method:
##$ -S /bin/sh
##$ -N assoc_test
##$ -cwd
##$ -j y
##$ -t 1-22
##Run with qsub -q new.q calculate_snp_pc_corr.sh

#module load conda

#R -q --vanilla --args CFF_sid_onlyGT.gds annot.rds 6-18pcr_mat.rds F508_count gaussian --out_prefix 6-18${SGE_TASK_ID} --covars "PC1 PC2 PC3" --variant_id var_filter_SNVs_MAF0.05.rds --sample_id keep_samples.rds --chromosome $SGE_TASK_ID < assoc_test.R &
