#!/bin/bash

for f in {1..22}; do
echo "Assoc_testing '$f' "
R -q --vanilla --args CFF_sid_onlyGT.gds annot.rds 6_26pcr_mat.rds F508_count poisson --out_prefix "F508del_noPCs_chr_$f" --covars "site" --variant_id SNPS_bi_GATK_VQSR_MAF0.05_miss0.05.rds --sample_id keep_samples_noTwins.rds --chromosome "$f" < assoc_test.R > F508del_assoc_test_noPCs.log &
sleep 30
done

#Run with: sh F508del_assoc_test_noPCs.sh
